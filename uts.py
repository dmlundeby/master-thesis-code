########################################
# Utility file for master project
#
# By Daniel Martin Lundeby, June 2019
########################################


import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import re
import numpy as np
import scipy
import hyperspy.api as hs
import datetime as dt
import time
import pickle
import json
from hyperspy.learn.mva import LearningResults
import copy
import pint
ureg = pint.UnitRegistry()
ureg.setup_matplotlib()


cfg = {
    'counts_per_electron': None    # TODO Remove this
}


path_common = '..\\common\\'


def example_signal_map():
    return hs.load(r'C:\Users\danie\Docsys-PC\Master_Data\signals\SCN45-ARM\maps\SL4.hspy')


def example_signal_spot():
    return hs.load(r'C:\Users\danie\Docsys-PC\Master_Data\signals\examples\A1.emsa', signal_type='EDS_TEM')


def save_report(size, *name):
    path_report = r'C:\Users\danie\Google Drive\Docsys-Gdrive\Utdanning\NTNU\Emner\Semester 12\TFY4905 Masteroppgave\Report-Latex\figures' + '\\'
    fix_fig(size)
    if len(name) == 0 or name[0] == '':
        raise Exception('You need to specify a name')
    plt.tight_layout()
    from os.path import join
    plt.savefig(join(path_report, ''.join(map(str,name))))


def fix_fig(width_fraction):
    fig = plt.gcf()
    fig.set_dpi(300)
    actual_width = 4.9616  # Width in inches of B5 paper used in project report\n",
    factor = 2
    if isinstance(width_fraction, float) or isinstance(int, float):
        aspect_ratio = 3 / 4
        wx = factor * actual_width * width_fraction
        wy = aspect_ratio * factor * actual_width * width_fraction
    else:
        wfx, wfy = width_fraction
        wx = factor * actual_width * wfx
        wy = factor * actual_width * wfy
    fig.set_size_inches(wx, wy)
    plt.tight_layout()


def load_files(path, str_format=None, file_types=None, signal_type=None, sort_by_ts=True, includeonly=None):
    """

    :param signal_type:
    :param file_types:
    :param path:
    :param str_format: Format of the filename to be recognized.
    Substitute {d} for number, {dd} for number which also supports numbers like `2-3` and {s} for string.
    :return:
    """

    files = []

    for name in os.listdir(path):

        # Ignore folders
        spl = name.rfind('.')
        if spl is -1:
            continue

        # Ignore wrong file types
        pre, post = name[:spl], name[spl + 1:]
        if file_types and post not in file_types:
            continue

        # Check format
        if str_format is not None:
            s = "[a-zA-zæøåÆØÅ_\-]+"
            d = "[0-9]+"
            dd = "[0-9]+(-[0-9]+)?"
            form = "^" + str_format.replace("{d}", d).replace("{s}", s).replace("{dd}", dd) + "$"
            if not re.match(form, pre):
                continue

        # Ignore if specified by includeonly
        if includeonly is not None and pre not in includeonly:
            continue

        filename = path + name
        s = hs.load(filename, signal_type=signal_type)
        s.metadata.General.set_item('title', pre)
        files.append(s)

    if sort_by_ts:
        files.sort(key=lambda s: get_ts(s))

    return files


def _get_all_bc_image(f, beam_current=True):
    """

    :param f:
    :param beam_current: If True, calculate beam current. If False, calculate counts per second.
    :return:
    """

    container = f.original_metadata.ImageList.TagGroup0.ImageTags
    exp = container.DataBar.Exposure_Time_s
    if 'Acquisition_series' in container:
        base = container.Acquisition_series

        cd = base.Sum.as_dictionary()
        cdkeys = sorted(cd.keys(), key=lambda x: int(x.split('_')[1]))
        # Cast to float in case the data is stored erroneously as a string (should be stored as float)
        counts = [float(cd[key]) for key in cdkeys]
        currents = counts_to_currents(counts, exp, beam_current)

        if 'Time' in base:
            td = base.Time.as_dictionary()
            times = [get_dt_gatan(td[key], f) for key in sorted(td.keys(), key=lambda x: int(x.split('_')[1]))]
        else:
            times = [get_dt(f) + dt.timedelta(seconds=3 * i) for i in range(len(currents))]

    else:
        currents = counts_to_currents([np.sum(f.data/1)], exp, beam_current)
        # For some weird reason, I need to divide the array data by 1. It is probably an issue with data formats.

        times = [get_dt(f)]

    return times, currents


def counts_to_currents(counts, exp, beam_current):
    if beam_current:
        if cfg['counts_per_electron'] is None:
            raise Exception('CCD conversion factor (X-ray counts per electron) needs to be set.')
        factor = (1 / cfg['counts_per_electron']) * scipy.constants.e * 1e9  # Unit: nC per count
        currents = [c / exp * factor for c in counts]
    else:
        currents = [c / exp for c in counts]
    return currents


def _get_bc_image(f, beam_current=True):
    tlist, clist = _get_all_bc_image(f, beam_current)
    current = np.median(clist)
    t1, t2 = get_ts(tlist[0]), get_ts(tlist[-1])
    time = get_dt( int(t1 + (t2 - t1)/2) )
    return time, current


def get_all_bc(files, beam_current=True):
    all_dts = []
    all_currents = []
    if isinstance(files, hs.signals.BaseSignal):
        files = [files]
    for file in files:
        dt, c = _get_all_bc_image(file, beam_current)
        all_dts.extend(dt)
        all_currents.extend(c)
    return all_dts, all_currents


def get_bc(files, beam_current=True):
    dts = []
    currents = []
    if isinstance(files, hs.signals.BaseSignal):
        files = [files]
    for f in files:
        dt, c = _get_bc_image(f, beam_current)
        currents.append(c)
        dts.append(dt)
    return dts, currents


def set_beam_current(files, signals, strategy='interp_fill', condition=None):
    """

    :param files:
    :param signals:
    :param strategy: One of `interp_fill`, `interp_extrap`, `regr` or `mean`.
    :param condition: Function with a signal as input, which determines whether to use the signal for estimating beam currents. This allows for excluding irregular beam current values.
    """

    ff = list(filter(condition, files))
    dts, currents = get_bc(ff)
    tss = [get_ts(dt) for dt in dts]

    get_beam_current = get_bc_func(currents, tss, strategy)    # Function for setting beam current

    # Set beam current values
    for s in signals:
        current = get_beam_current(get_ts(s))  # Retrieve beam current value
        s.set_microscope_parameters(beam_current=current)  # in nA


def get_bc_func(currents, tss, strategy):
    # Interpolate
    if strategy == 'interp_extrap':
        c_int = scipy.interpolate.interp1d(tss, currents, kind='linear', fill_value='extrapolate')  # Interpolated curve
        return lambda x: c_int(x).item()

    if strategy == 'interp_fill':
        int_func = scipy.interpolate.interp1d(tss, currents, kind='linear')  # Interpolated curve
        def c_int(x):  # Fills values
            try:
                return int_func(x).item()
            except:
                if x <= tss[0]:
                    return currents[0]
                else:
                    return currents[-1]
        return c_int

    if strategy == 'regr':
        regr = scipy.stats.linregress(tss, currents)
        return lambda x: regr.intercept + x * regr.slope

    if strategy == 'mean':
        return lambda x: np.mean(currents)

    if strategy == 'middle':
        def middle(x):
            index = -1
            for i, v in enumerate(tss):
                if x >= v:
                    index = i
            if index == -1:
                return currents[0]
            elif index == len(tss) - 1:
                return currents[-1]
            else:
                return (currents[index] + currents[index + 1]) / 2
        return middle

    raise Exception('The strategy `' + strategy + '` does not exist.')

def gets(signals, name):
    return next(filter(lambda s: s.metadata.General.title == name, signals))


def get_ts(s):
    """

    :param s: Either a BaseSignal or a datetime object.
    :return:
    """
    if isinstance(s, hs.signals.BaseSignal):
        s = get_dt(s)
    if isinstance(s, dt.datetime):
        return int(time.mktime(s.timetuple()))


def get_dt(s):
    """

    :param s: Either a BaseSignal or a timestamp integer.
    :return:
    """
    if isinstance(s, hs.signals.BaseSignal):
        date = s.metadata.General.date
        time = s.metadata.General.time
        datetime = dt.datetime.strptime(date + ' ' + time, '%Y-%m-%d %H:%M:%S')
        return datetime
    elif isinstance(s, int):
        return dt.datetime.fromtimestamp(s)


def get_dt_gatan(datetime_str, file=None):
    # Should be implemented better
    fmt = '%d/%m/%Y %I:%M:%S %p'
    if file:
        microscope_info = file.original_metadata.ImageList.TagGroup0.ImageTags.Microscope_Info
        if 'Items' in microscope_info:
            microscope = microscope_info.Items.TagGroup2.Value
            if microscope == '2100':
                fmt = '%m/%d/%Y %I:%M:%S %p'
    return dt.datetime.strptime(datetime_str, fmt)


def _split(text):
    if isinstance(text, hs.signals.BaseSignal):
        text = get_label(text)

    index = 0
    for i, t in enumerate(text):
        if t.isdigit():
            index = i
            break
    return text[:index], text[index:]


def get_label(s):
    return s.metadata.General.label

def get_region(s):
    return _split(s)[0]


def get_index(s):
    return _split(s)[1]


def get_first_index(s):
    return get_index(s).split('-')[0]


def create_tilt_series(signals, regions, xray_lines, get_intensities, get_composition, get_thickness, get_density, exclude=[]):
    tilt_series = []

    def tiltx(s):
        return s.metadata.Acquisition_instrument.TEM.Stage.tilt_alpha

    for line in xray_lines:
        for reg in regions:
            ss = [s for s in signals if s.metadata.General.title.startswith(reg)]

            def filter_exclude(x):
                label, region, index, first_index = get_label(x), get_region(x), get_index(x), get_first_index(x)
                for exc in exclude:
                    if exc == label or exc == region or exc == index or exc == first_index:
                        return False
                return True
            ss = list(filter(filter_exclude, ss))

            ss.sort(key=lambda s: tiltx(s))
            reg_name = ss[0].metadata.General.region_name
            els = get_composition(ss[0]).keys()
            tilts = np.array([tiltx(s) for s in ss])
            labels = np.array([s.metadata.General.title for s in ss])

            el = line.split('_')[0]
            if el in els:
                zetas = []
                for s in ss:
                    zfac = determine_zeta_factor(
                        s,
                        [line],
                        get_intensities,
                        get_composition,
                        get_thickness,
                        get_density)[0].data[0]
                    zetas.append(zfac)

                tilt_series.append({
                    'tilts': np.array(tilts),
                    'zetas': np.array(zetas),
                    'xray_line': line,
                    'labels': labels,
                    'region': reg,
                    'region_name': reg_name})

    return tilt_series


def _get_color(prop, colors):
    color_list = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
    if not prop in colors:
        for c in color_list:
            if c not in colors.values():
                colors[prop] = c
                break
    return colors[prop]


def _get_marker(line, marker):
    markers = ['*', '.']
    el = line.split('_')[0]
    if not el in marker:
        for m in markers:
            if m not in marker.values():
                marker[el] = m
    return marker[el]


def plot_tilt_series(tilt_series, regions=None, xray_lines=None, avg=False,
                     annotate=False, ann_id=0, ann_off=None, lines=True, zone_axes=False, zone_circles=False,
                     plot_exp_fit=False,
                     xmin=None, xmax=None, ymin=None, ymax=None, ref=None,
                     fig=None, ax=None, legend=True, color='multi', ax_deg=True):
    if fig is None and ax is None:
        fig, ax = plt.subplots()
    col = {}
    marker = {}

    tilt_series = ts_filter(tilt_series, regions, xray_lines, avg)

    for i, ts in enumerate(tilt_series):
        if legend is True:
            L = format_legend(ts, '{l}{,rn}')
        elif type(legend)==str:
            L = format_legend(ts, legend)
        elif type(legend)==list:
            L = format_legend(ts, legend[i])
        else:
            L = None
        m = _get_marker(ts['xray_line'], marker)
        if color == 'multi':
            c = _get_color(ts['region'], col)
        elif color == 'iter':
            c = None
        else:
            c = color
        linestyle = '-' if lines else ''
        ax.plot('tilts', 'zetas', data=ts, marker=m, color=c, linestyle=linestyle, label=L)

        if plot_exp_fit:
            res, cov = fit_shadowing(ts)
            plot_shadowing(res, ts, ax=ax)


    if ymin is None: ymin = 0
    if xmin is not None: plt.xlim(left=xmin)
    if xmax is not None: plt.xlim(right=xmax)
    if ymin is not None: plt.ylim(bottom=ymin)
    if ymax is not None: plt.ylim(top=ymax)

    diff = ax.get_ylim()[1] - ax.get_ylim()[0]
    if annotate not in [None, False]:
        if annotate=='all':
            if ann_off == None: ann_off = 0
            for ts in tilt_series:
                for tilt, zeta, label in zip(ts['tilts'], ts['zetas'], ts['labels']):
                    if in_range(tilt, zeta+ann_off*diff, ax):
                        ax.annotate(label, (tilt, zeta+ann_off*diff))
        elif annotate=='single':
            if ann_off == None: ann_off = -0.1
            ts = tilt_series[ann_id]
            for tilt, zeta, label in zip(ts['tilts'], ts['zetas'], ts['labels']):
                if in_range(tilt, zeta+ann_off*diff, ax):
                    ax.annotate(get_first_index(label), (tilt, zeta+ann_off*diff))
        else:
            raise Exception('Unknown annotation option: ' + str(annotate))

    if zone_axes is not False:
        plot_zone_axes(tilt_series[0]['zone_axes'], annotate=True)

    if zone_circles is not False:
        for t in tilt_series:
            for i, tilt in enumerate(t['tilts']):
                if tilt in t['zone_axes']:
                    plot_circle(tilt, t['zetas'][i], t['zone_axes'][tilt])

    if legend: ax.legend()

    _plot_tilt_series(fig, ax, ax_deg)
    if ref is not None:
        add_deviation(ax, ref)
    return fig, ax


def format_legend(ts, text):
    if text is None:
        return None
    l = ts['xray_line']
    e = l.split('_')[0]
    r = ts['region']
    rn = ts['region_name']
    inst = ts['instrument'] if 'instrument' in ts else ''
    set_id = ts['id'] if 'id' in ts else ''
    text = leg_repl(text, '{e}', e)
    text = leg_repl(text, '{l}', _style_line(l))
    text = leg_repl(text, '{r}', r)
    text = leg_repl(text, '{rn}', rn)
    text = leg_repl(text, '{i}', inst)
    text = leg_repl(text, '{id}', set_id)
    return text


def leg_repl(text, repl, item):
    text = text.replace(repl, item)

    repl_comma = repl[0]+','+repl[1:]
    if item is not None and item is not '':
        item = ', ' + item
    return text.replace(repl_comma, item)


def ts_filter(ts, regions=None, xray_lines=None, avg=None, exclude=None):
    ts = ts_reduce(ts, regions, xray_lines)
    ts = ts_avg(ts, avg)
    ts = ts_exclude(ts, exclude)
    return ts


def ts_reduce(ts, regions=None, xray_lines=None):
    ts_out = []
    for t in ts:
        if (
                (regions is None or t['region'] in regions) and
                (xray_lines is None or t['xray_line'] in xray_lines)
        ):
            ts_out.append(t)
    return ts_out


def ts_avg(ts, option='xray_lines'):

    if option == 'all':
        return _ts_avg(ts)

    elif option == 'xray_lines':
        ts_out = []
        xray_lines = ts_get_xray_lines(ts)
        for line in xray_lines:
            tsr = ts_reduce(ts, xray_lines=[line])
            ts_out.extend(_ts_avg(tsr))
        return ts_out

    elif option is False or option is None:
        return ts

    else:
        raise Exception('Unknown averaging option: ' + str(option) + '. Allowed options are [`all`, `xray_lines`]')


def ts_exclude(ts, indices=None):
    if indices is None: return ts
    ts_out = copy.deepcopy(ts)
    for t in ts_out:
        tilts = []
        zetas = []
        labels = []
        for tilt, zeta, label in zip(t['tilts'], t['zetas'], t['labels']):
            index, first_index = get_index(label), get_first_index(label)
            if index not in indices and first_index not in indices:
                tilts.append(tilt)
                zetas.append(zeta)
                labels.append(label)
        t['tilts'] = np.array(tilts)
        t['zetas'] = np.array(zetas)
        t['labels'] = np.array(labels)
    return ts_out

def ts_get_xray_lines(ts):
    xl = set()
    for t in ts:
        xl.add(t['xray_line'])
    return list(xl)


def _ts_avg(ts):
    zetas = np.zeros(len(ts[0]['zetas']))
    for t in ts:
        zetas += t['zetas']
    zetas /= len(ts)
    ts_dict = {'tilts': ts[0]['tilts'], 'zetas': zetas, 'xray_line': ts[0]['xray_line'], 'region': '',
            'region_name': '', 'labels': ts[0]['labels'], 'zone_tilt': ts[0]['zone_tilt'],
             'zone_axes': ts[0]['zone_axes']}
    if 'id' in ts[0]:
        ts_dict.update({'id': ts[0]['id']})
    return [ts_dict]


def in_range(x, y, ax):
    xlim = ax.get_xlim()
    if x < xlim[0] or x > xlim[1]:
        return False
    ylim = ax.get_ylim()
    if y < ylim[0] or y > ylim[1]:
        return False
    return True


def plot_zone_tilts(ax, zone_tilt):
    hits = set()
    tmp = zone_tilt
    while tmp >= ax.get_xlim()[0]:
        hits.add(tmp)
        tmp = tmp - 30
    tmp = zone_tilt
    while tmp <= ax.get_xlim()[1]:
        hits.add(tmp)
        tmp = tmp + 30
    for hit in hits:
        plt.axvline(x=hit, linestyle='--', color='grey')


def plot_zone_axes(zone_axes, annotate=False):
    for tilt, text in zone_axes.items():
        plt.axvline(x=tilt, linestyle='--', color='grey')

        if annotate:
            ymin, ymax = plt.ylim()
            ydiff = ymax-ymin
            xmin, xmax = plt.xlim()
            xdiff = xmax-xmin
            plt.annotate('[' + text + ']', (tilt + 0.01 * xdiff, ymin + 0.9 * ydiff), rotation=90)


def plot_circle(x, y, text=None, xfac=0.02, yfac=0.01):
    ymin, ymax = plt.ylim()
    ydiff = ymax - ymin
    xmin, xmax = plt.xlim()
    xdiff = xmax - xmin

    if in_range(x, y, plt.gca()):
        plt.scatter([x], [y], s=100, edgecolor='k', facecolor='none', zorder=10)
        plt.annotate('[' + text + ']', (x + xfac*xdiff, y - yfac*ydiff))


def circle_zeta(ts, text, label, xray_line=None, regions=None, xray_lines=None, avg=False):

    ts = ts_reduce(ts, regions, xray_lines)
    ts = ts_avg(ts, avg)

    i = ts_getzeta(ts, label, xray_line)
    if i:
        plot_circle(i['tilt'], i['zeta'], text)


def ts_getzeta(ts, label, xray_line=None):
    for t in ts:
        if xray_line is None or t['xray_line'] == xray_line:
            if label in t['labels']:
                i = int(np.where(t['labels']==label)[0][0])
                return {'tilt': t['tilts'][i], 'zeta': t['zetas'][i]}


def to_hours(dts, ref=None):
    if ref is None:
        ref = dts[0] if isinstance(dts, list) else dts
    f = lambda t: ((t - ref).total_seconds())/3600
    if isinstance(dts, list):
        return [f(t) for t in dts]*ureg.hour
    else:
        return f(dts)*ureg.hour


def to_minutes(dts, ref=None):
    if ref is None:
        ref = dts[0] if isinstance(dts, list) else dts
    f = lambda t: ((t - ref).total_seconds())/60
    if isinstance(dts, list):
        return [f(t) for t in dts]*ureg.minute
    else:
        return f(dts)*ureg.minute


def to_rel_time(dts, ref=None):
    plt.gca().xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))
    if ref is None:
        ref = dts[0] if isinstance(dts, list) else dts
    f = lambda t: dt.datetime.fromtimestamp(0) - dt.timedelta(hours=1) + (t - ref)
    if isinstance(dts, list):
        return [f(t) for t in dts]
    else:
        return f(dts)


def to_abs_time(dts, ref=None):
    plt.gca().xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))
    return dts


def _plot_tilt_series(fig, ax, ax_deg=True):
    if ax_deg:
        ax.xaxis.set_major_formatter(mpl.ticker.StrMethodFormatter(u'{x:.0f}°'))
        ax.set_xlabel('Stage tilt')
    else:
        ax.set_xlabel('Stage tilt (degrees)')
    ax.set_ylabel(r'$\zeta$-factor (kg electron m$^{-2}$ photon$^{-1}$)')
    fig.tight_layout()


def plot_beam_current(files, mean=True, full=False, signals=None, conv=to_hours, ref=None, ref_dec=0, annotate='single',
                      xs='time', label=True, beam_current=True, linestyle='-'):
    dts, cs = get_bc(files, beam_current)
    adts, acs = get_all_bc(files, beam_current)

    if beam_current is False:
        cs = [c/1e9 for c in cs]
        acs = [c/1e9 for c in acs]

    fig, ax = plt.subplots()

    reftime = dts[0]
    if signals is not None:
        signal_min = min([get_dt(s) for s in signals])
        if signal_min < dts[0]:
            reftime = signal_min

    # Plot beam current measurements
    if mean:
        if xs == 'time':
            y = conv(dts, ref=reftime)
        else:
            if xs not in files[0].metadata: raise Exception('Entry ' + xs + ' is not in metadata')
            y = [f.metadata.get_item(xs) for f in files]
        ax.plot(y, cs, marker='o', linestyle=linestyle, label='Measured' if label else None)

    # Plot all beam current data
    if full:
        if xs == 'time':
            y = conv(adts, ref=reftime)
        else:
            if xs not in files[0].metadata: raise Exception('Entry ' + xs + ' is not in metadata')
            y = []
            for f in files:
                y.extend([f.metadata.get_item(xs)] * len(f.original_metadata.ImageList.TagGroup0.ImageTags.Acquisition_series.Time))
        ax.scatter(y, acs, c='g', s=10, label='Individual measurements' if label else None)

    # Plot interpolated beam currents, for each EDS spectrum
    if signals:
        x = [conv([get_dt(s)], ref=reftime) for s in signals]
        y = [s.metadata.Acquisition_instrument.TEM.beam_current for s in signals]
        ax.scatter(x, y, c='red', s=20, label='Interpolated' if label else None)

        if annotate not in [None, False]:
            if annotate=='all':
                sgns = signals
            elif annotate=='single':
                sgns = get_firstindex_signals(signals)
            else:
                raise Exception('Unknown annotation option: ' + str(annotate))
            for s in sgns:
                xi = conv([get_dt(s)], ref=reftime)
                yi = s.metadata.Acquisition_instrument.TEM.beam_current
                if annotate=='all': text = get_label(s)
                elif annotate=='single': text = get_first_index(s)
                ax.annotate(text, (xi, yi))

    if ref is not None:
        add_deviation(ax, ref, ref_dec)

    #fig.autofmt_xdate()
    if xs == 'time': ax.set_xlabel('Time elapsed (' + ax.get_xlabel() + 's)')
    if beam_current: ax.set_ylabel('Beam current (nA)')
    else: ax.set_ylabel('Count rate ($10^9$ cps)')
    if label: ax.legend()
    fig.tight_layout()
    return fig, ax


def get_firstindex_signals(signals):
    out = []
    indices = set()
    for s in signals:
        index = get_first_index(s)
        if index not in indices:
            out.append(s)
        indices.add(index)
    return out


def ts_get(ts, region, xray_line):
    return next(filter(lambda x: x['region'] == region and x['xray_line'] == xray_line, ts))


def filter_zetas(ts, indices):
    if indices is None:
        return ts['zetas']
    if isinstance(indices, dict):
        indices = indices[ts['id']]
    out = []
    for label in ts['labels']:
        if get_index(label) in indices:
            out.append(True)
        else:
            out.append(False)
    return ts['zetas'][out]


def get_zeta_factors(tilt_series, indices=None, regions=None, xray_lines=None, avg=False):

    tilt_series = ts_filter(tilt_series, regions, xray_lines, avg)

    zeta_factors = []

    for ts in tilt_series:

        zetas = filter_zetas(ts, indices)

        zdict = {'xray_line': ts['xray_line'], 'region': ts['region'], 'region_name': ts['region_name'],
             'zeta': np.mean(zetas), 'std': np.std(zetas)}
        if 'id' in ts:
            zdict.update({'id': ts['id']})
        zeta_factors.append(zdict)

    # Might add an option for adding "total zeta-factors" for a given element.

    return zeta_factors


def ts_merge(ts, indices=None, regions=None, xray_lines=None, avg=False):

    ts = ts_filter(ts, regions, xray_lines, avg)

    ts0 = copy.deepcopy(ts[0])
    zs = []
    for t in ts:
        zs.extend(list(filter_zetas(t, indices)))
        # TODO extend remaining variables
    ts0['zetas'] = np.array(zs)
    return [ts0]




def plot_zeta_factors(zeta_factors, color='region_name', label='{l}{,rn}{,id}', xfac=0):

    fig, ax = plt.subplots()
    colors = {}
    marker = {}

    for zf in zeta_factors:
        L = format_legend(zf, label)
        m = _get_marker(zf['xray_line'], marker)
        if color == 'iter':
            c = None
        elif color in zf:
            c = _get_color(zf[color], colors)
        else:
            c = color

        ax.errorbar([L], zf['zeta'], zf['std'], marker=m, color=c, linestyle='', capsize=3, ms=10)

    _plot_zeta_factors(fig, ax, xfac)
    return fig, ax


def _style_line(line):
    ind_bar = line.index('_')
    outline = line[:ind_bar] + ' ' + line[ind_bar+1] + r'$\mathregular{_'
    out = None
    lt = line[ind_bar+2]
    if lt == 'a': out = u'α'
    elif lt == 'b': out = u'β'
    elif lt == 'g': out = u'γ'
    elif lt == 'n': out = u'μ'
    elif lt == 'l': out = u'ℓ'
    outline += out + r'}$'
    if len(line)>ind_bar+3:
        outline += line[ind_bar+3:]
    return outline

def _plot_zeta_factors(fig, ax, xfac=0):
    ax.set_ylabel(r'Zeta-factor (kg electron m$^{-2}$ photon$^{-1}$)')

    lm = ax.get_xlim()
    mn, mx = lm[0], lm[1]
    diff = mx - mn
    ax.set_xlim(mn - xfac * diff, mx + xfac * diff)

    plt.xticks(rotation=45)

    add_deviation(ax)

    fig.tight_layout()


def plot_tilt_series_comparison(tss, labels, region, xray_lines):
    fig, ax = plt.subplots()
    col = {}
    marker = {}
    for ts, label in zip(tss, labels):
        for line in xray_lines:
            m = _get_marker(line, marker)
            c = _get_color(label, col)
            data = ts_get(ts.tilt_series, region, line)
            ax.plot('tilts', 'zetas', data=data, marker=m, color=c, label=line + ', ' + label)
    ax.legend()
    _plot_tilt_series(fig, ax)
    return fig, ax


def add_deviation(ax, ref=None, ref_dec=0):
    ax2 = ax.twinx()
    mn, mx = ax.get_ylim()
    if ref is None:
        ref = (mn + mx)/2
    ax2.set_ylim((mn - ref) / ref, (mx - ref) / ref)
    ax2.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:.' + str(ref_dec) + '%}'))
    ax2.set_ylabel('Deviation')
    plt.gcf().tight_layout()


def _get_param_str(p, val):
    if p.name == 'A':
        if val is not None:
            return '{:<10.0f}'.format(val)
        else:
            return '{:10s}'.format('')
    elif p.name == 'centre':
        if val is not None:
            return '{:<10.3f}'.format(val)
        else:
            return '{:10s}'.format('')
    elif p.name == 'sigma':
        if val is not None:
            return '{:<10.3g}'.format(val)
        else:
            return '{:10s}'.format('')


def print_model(m):
    s = '{:11s}'.format('Component')
    s += '{:10s}'.format('Name')
    s += '{:10s}'.format('Value')
    s += '{:10s}'.format('Min')
    s += '{:10s}'.format('Max')
    s += '{:10s}'.format('Twin')
    s += '{:15s}'.format('Free')
    s += '\n'

    for c in m:
        if c.name.startswith('background'): continue
        s += c.name + '\n'
        for p in c.parameters:
            s += '{:11s}'.format('')
            s += '{:10s}'.format(p.name)
            s += _get_param_str(p, p.value)
            s += _get_param_str(p, p.bmin)
            s += _get_param_str(p, p.bmax)
            s += '{:10s}'.format(p.twin.component.name + ', ' + p.twin.name[0] if p.twin else '')
            s += '{:15s}'.format('T' if p.free else 'F')
            s += '\n'
        s += '\n'
    print(s)


def _init_window_var(xray_lines, var, default):
    if isinstance(var, dict): return var
    val = var if var is not None else default
    return {line: val for line in xray_lines}


def get_bws(s, xray_lines, bw_pos, bw_width):
    bw_pos = _init_window_var(xray_lines, bw_pos, [2, 2])
    bw_width = _init_window_var(xray_lines, bw_width, 1)
    bw = []
    for line in xray_lines:
        bw.append(
            s.estimate_background_windows(xray_lines=[line], windows_width=bw_width[line], line_width=bw_pos[line])[0]
        )
    return np.array(bw)


def get_iws(s, xray_lines, iw_width):
    iw_width = _init_window_var(xray_lines, iw_width, 2)
    iw = []
    for line in xray_lines:
        iw.append(
            s.estimate_integration_windows(xray_lines=[line], windows_width=iw_width[line])[0]
        )
    return iw    # Causes a bug if made as numpy array


def get_int_win(s, xray_lines, bw_pos, bw_width, iw_width):
    # xray_lines are really contained within bw_pos
    bw = get_bws(s, xray_lines, bw_pos, bw_width)
    iw = get_iws(s, xray_lines, iw_width)
    li = s.get_lines_intensity(xray_lines=xray_lines, integration_windows=iw, background_windows=bw)
    return li


def plot_window_results(s, xray_lines, bw_pos, bw_width, iw_width):
    bw = get_bws(s, xray_lines, bw_pos, bw_width)
    iw = get_iws(s, xray_lines, iw_width)
    s.plot(xray_lines='from_elements', background_windows=bw, integration_windows=iw)


def _extract_info(elem):

    if callable(elem):
        return elem

    if not isinstance(elem, dict):
        return lambda x: elem

    def _get_value(item, s):
        if callable(item):
            return item(s)
        else:
            return item

    def _extract_from_dict(s):

        # Priority order: Label --> Region --> Index

        label = get_label(s)
        if label in elem:
            return _get_value(elem[label], s)
        region = get_region(s)
        if region in elem:
            return _get_value(elem[region], s)
        index = get_index(s)
        if index in elem:
            return _get_value(elem[index], s)
        else:
            raise Exception('Cannot find any matching element for the following signal: ' + get_label(s) +
                            ', in the element ' + str(elem))

    return _extract_from_dict


def determine_zeta_factor(s, xray_lines, get_intensities, get_composition, get_thickness, get_density, get_dose=None):
    if isinstance(get_intensities, dict):
        intensities = [get_intensities[get_label(s)][line][0] for line in xray_lines]
    else:    # Function
        intensities = get_intensities(s, xray_lines)
    composition = [get_composition(s)[line.split('_')[0]] for line in xray_lines]
    thickness = get_thickness(s)
    density = get_density(s)
    if get_dose is None:
        dose = s._get_dose('zeta')
    else:
        dose = get_dose(s)

    return _determine_zeta_factor(intensities, composition, thickness, density, dose)


def determine_relative_zeta_factor(s, xray_lines, get_int, get_comp):
    return determine_zeta_factor(s, xray_lines, get_int, get_comp, lambda x: 1, lambda x: 1, lambda x: 1)


def determine_zeta_factor_from_internal(s, xray_lines, get_internal_zeta, get_thickness, get_density, get_dose=None):
    comp_construct = lambda x: {xray_line.split('_')[0]: get_internal_zeta(x) for xray_line in xray_lines}
    int_construct = lambda x, y: [1 for xray_line in xray_lines]
    return determine_zeta_factor(s, xray_lines, int_construct, comp_construct, get_thickness, get_density, get_dose)


def _determine_zeta_factor(intensities, composition, thickness, density, dose):
    return [np.divide(dose*density*thickness*comp, intensity) for comp, intensity in zip(composition, intensities)]


class BeamCurrentAnalyzer:

    def __init__(self, params):
        self.params = params
        self.folder_current = params['folder_current'].rstrip('\\') + '\\'
        self.load_beam_currents(params)

    def load_beam_currents(self, params):
        cfg['counts_per_electron'] = params['counts_per_electron']  # Should fix this
        form = self.params['format'] if 'format' in self.params else 'current{dd}'
        sbt = params['sort_by_time'] if 'sort_by_time' in self.params else True

        inconly = self.params['includeonly'] if 'includeonly' in params else None
        self.beam_currents = load_files(self.folder_current, form, file_types=['dm3', 'dm4'], sort_by_ts=sbt,
                                        includeonly=inconly)


    def get_beam_current(self, beam_current=True):
        if hasattr(self, 'beam_currents'):
            return get_bc(self.beam_currents, beam_current)
        else:
            return self.params['beam_current']

    def plot_beam_current(self, mean=True, full=False, interpolated=True, conv=to_hours, ref=None, ref_dec=0,
                          annotate='single', xs='time', label=True, beam_current=True):
        if hasattr(self, 'beam_currents'):
            sgns = self.signals if (interpolated and hasattr(self, 'signals')) else None
            return plot_beam_current(self.beam_currents, mean, full, sgns, conv, ref, ref_dec, annotate, xs, label,
                                     beam_current)
        else:
            fig, ax = plt.subplots()
            ax.plot(self.params['beam_current'].keys(), self.params['beam_current'].values(), 'o')
            return ax, fig


class ModelFitter:

    signals = []
    models = []
    models_sum = []
    get_int = []

    def __init__(self, signals, label=None):
        self.signals = signals
        self.path_model_cache = '.\\cache_model' + (('_' + label) if label else '') + '\\'
        self.path_model_cache_models = self.path_model_cache + 'models\\'
        self.path_model_cache_models_sum = self.path_model_cache + 'models_sum\\'
        self.path_model_intensities = self.path_model_cache + 'intensities.json'

    def model_fit(self, mod_lim, mod_elements=None, mod_xray_lines=None, use_cache=True, load_models=False):

        self.models = []
        self.models_sum = []
        self.get_int = {}

        if use_cache and self._has_model_cache():
            print('Loading model fitting results from cache')
            if load_models:
                self._load_model_cache()
            try:
                self._load_ints_from_cache()
            except Exception:
                print('Error loading intensities from cache')
        else:
            for s in self.signals:
                print('Model fitting signal ' + s.metadata.General.title)
                s, m, m_sum = self._model_fit(s, mod_lim, mod_elements, mod_xray_lines)
                self.models.append(m)
                self.models_sum.append(m_sum)
                self._create_model_cache(m, m_sum)
            self._create_intensity_cache()
            try:
                self._load_ints_from_cache()
            except Exception:
                print('Error loading intensities from cache')

    def model_fit_sum(self, mod_lim, mod_elements=None, mod_xray_lines=None, label=None):
        s = self._mod_get_s(label).sum()
        s, m, m_sum = self._model_fit(s, mod_lim, mod_elements, mod_xray_lines)
        m.plot(plot_components=True, xray_lines='from_elements', only_lines=None)

    def plot_model(self, label=None):
        m = self._mod_get_m(label)
        m.plot(plot_components=True, xray_lines='from_elements', only_lines=None)

    def _create_model_cache(self, m, m_sum):
        m.signal.learning_results = LearningResults() # Avoid HyperSpy error

        if not os.path.isdir(self.path_model_cache_models):
            os.makedirs(self.path_model_cache_models)
        m.save(self.path_model_cache_models + m.signal.metadata.General.title + '.hspy', overwrite=True)

        if m_sum is not None:
            m_sum.signal.learning_results = LearningResults() # Avoid HyperSpy error
            if not os.path.isdir(self.path_model_cache_models_sum):
                os.makedirs(self.path_model_cache_models_sum)
            m_sum.save(self.path_model_cache_models_sum + m.signal.metadata.General.title + '.hspy', overwrite=True)

    def _model_fit(self, s, mod_lim, mod_elements=None, mod_xray_lines=None, **kwargs):
        multidim = len(s.axes_manager.navigation_axes) != 0

        s = s.isig[_extract_info(mod_lim)(s)]  # Limit the signal range
        s.metadata.set_item('Sample.elements', [])
        s.metadata.set_item('Sample.xray_lines', [])
        if mod_xray_lines:
            s.set_lines(_extract_info(mod_xray_lines)(s))
        elif mod_elements:
            s.set_elements(_extract_info(mod_elements)(s))
        else:
            raise Exception('Either `mod_xray_lines` or `mod_elements´ needs to be set')

        m_sum = None

        if multidim:
            m_sum = s.sum().create_model()
            self._remove_twins(m_sum)
            m_sum.set_parameters_free()
            self._set_bounds(m_sum)
            m_sum.fit(bounded=True)
            mp = self._model_params(m_sum)
            m = s.create_model()
            self._load_param_from_dict(m, mp)
            m.multifit(bounded=True, max_nfev=50, **kwargs)
            # max_nfev: Limits the amount of time to do model fitting, to avoid using an unreasonable amount of time.
            # (it seems that the model fittings gets stuck at some point in time)
        else:
            m = s.create_model()
            self._remove_twins(m)  # Remove all twins
            m.set_parameters_free()  # Free all parameters
            self._set_bounds(m)  # Set bounds as defined by this method
            m.fit(bounded=True, max_nfev=50, **kwargs)

        return s, m, m_sum

    def _remove_twins(self, m):
        for c in m:
            if c.name.startswith('background'):
                continue
            c.A.twin = None

    def _set_bounds(self, m, centre_offset=0.01, sigma_offset=0.1):
        for c in m:
            if c.name.startswith('background'): continue

            A = c.parameters[0]
            A.bmin = 0

            centre = c.parameters[1]
            centre.bmin = centre.value - centre.value * centre_offset
            centre.bmax = centre.value + centre.value * centre_offset

            sigma = c.parameters[2]
            sigma.bmin = sigma.value - sigma.value * sigma_offset
            sigma.bmax = sigma.value + sigma.value * sigma_offset

    def _load_model_cache(self):

        # Load models
        for fn in os.listdir(self.path_model_cache_models):
            m = hs.load(self.path_model_cache_models + fn)
            self.models.append(m.models.restore('a'))

        # Load model sums
        if os.path.isdir(self.path_model_cache_models_sum):
            for fn in os.listdir(self.path_model_cache_models_sum):
                m_sum = hs.load(self.path_model_cache_models_sum + fn)
                self.models_sum.append(m_sum.models.restore('a'))

    def _create_intensity_cache(self):
        out = {}
        for m in self.models:
            title = m.signal.metadata.General.title
            ints = {i.metadata.Sample.xray_lines[0]: [i] for i in m.get_lines_intensity()}
            out[title] = ints
        out_save = copy.deepcopy(out)
        for label, v in out_save.items():
            for xray_line, arr in v.items():
                if arr is not None:
                    out_save[label][xray_line] = [np.ndarray.tolist(s.data) for s in arr]
        with open(self.path_model_intensities, 'w') as file:
            json.dump(out_save, file)

    def _load_ints_from_cache(self):
        with open(self.path_model_intensities, 'r') as file:
            d = json.load(file)
        inp = copy.deepcopy(d)
        for label, v in d.items():
            signal = self._mod_get_s(label)
            for xray_line, arr in v.items():
                if arr is not None:
                    if isinstance(arr, list): arr = arr[0]
                    img = hs.signals.BaseSignal(signal.deepcopy())
                    img.metadata = signal.metadata.deepcopy()
                    img.axes_manager = signal.axes_manager.deepcopy()
                    am = img.axes_manager
                    am.remove(am.signal_axes)
                    if len(am.navigation_axes) == 0:
                        img.axes_manager.create_axes([{'size': 1, 'name': 'Scalar'}])
                    img.metadata.General.title = 'Intensity of ' + xray_line
                    img.metadata.set_item('Sample.elements', [xray_line.split('_')[0]])
                    img.metadata.set_item('Sample.xray_lines', [xray_line])
                    img.data = np.array(arr)
                    inp[label][xray_line] = [img]
        self.get_int = inp

    def _mod_get_s(self, label):
        if label is None:
            return self.signals[0]
        for s in self.signals:
            if s.metadata.General.title == label:
                return s

    def _mod_get_m(self, label):
        if label is None:
            return self.models[0]
        for m in self.models:
            if m.signal.metadata.General.title == label:
                return m

    def _has_model_cache(self):
        return os.path.isdir(self.path_model_cache)

    def _model_params(self, m):
        p_out = {}
        for c in m:
            if c.name.startswith('background'): continue
            name = c.name
            centre = c.centre.value
            sigma = c.sigma.value

            twin = getattr(m.components, name[:4] + 'a')
            twin_name = twin.name
            if name != twin_name:
                A_twin_ratio = c.A.value / twin.A.value
            else:
                A_twin_ratio = None
                twin_name = None
            p_out[name] = {'centre': centre, 'sigma': sigma, 'A_twin': twin_name, 'A_twin_ratio': A_twin_ratio}
        return p_out

    def _load_param_from_dict(self, m, d):
        for line, val in d.items():
            m.set_parameters_value('centre', val['centre'], [line])
            m.set_parameters_value('sigma', val['sigma'], [line])

            if val['A_twin'] is not None:
                c = getattr(m.components, line)
                twin = getattr(m.components, val['A_twin'])
                c.A.twin = twin.A
                c.A.twin_function_expr = 'x * ' + str(val['A_twin_ratio'])


class Analyzer(BeamCurrentAnalyzer, ModelFitter):

    def __init__(self, params, calibration=None):
        self.initialize(params, calibration)

    def initialize(self, params, calibration=None):
        self.params = params
        self.folder = params['folder']
        if not self.folder.endswith('\\'):
            self.folder = self.folder + '\\'
        self.folder_current = self.folder + 'beam_current\\'
        self._load_signals(params)

        for s in self.signals:
            if 'spectrum_calibration_scale' in params:
                params['spectrum_calibration_scale'] = s.axes_manager.signal_axes[0].scale
            if 'spectrum_calibration_offset' in params:
                params['spectrum_calibration_offset'] = s.axes_manager.signal_axes[0].offset
            s.set_elements(params['all_elements'])
            self.xray_lines = params['xray_lines']
            if calibration is not None:
                calibration(s)

        if 'beam_current' in params:
            if params['beam_current'] == 'from_measurement_files':
                self.load_beam_currents(params)
                strategy = params['beam_current_strategy'] if 'beam_current_strategy' in params else 'interp_fill'
                condition = params['beam_current_condition'] if 'beam_current_condition' in params else None
                set_beam_current(self.beam_currents, self.signals, strategy=strategy, condition=condition)

            else:
                get_beam_current = _extract_info(params['beam_current'])
                for s in self.signals:
                    s.set_microscope_parameters(beam_current=get_beam_current(s))  # in nA

        self.bw_pos = self._ext('bw_pos')
        self.bw_width = self._ext('bw_width')
        self.iw_width = self._ext('iw_width')
        if params['intensity'] == 'window_method':
            self._do_window_method()

        elif self.params['intensity'] == 'model_fitting':
            ModelFitter.__init__(self, self.signals, self.params['id'])
            self.model_fit(params['mod_lim'], self._ext('mod_elements'), self._ext('mod_xray_lines'),
                           use_cache=params['mod_cache'] if 'mod_cache' in params else True,
                           load_models=params['load_models'] if 'load_models' in params else False)

        def _comp_weight_percent(x, comp_atomic_percent):
            comp = comp_atomic_percent(x)
            weight_percent = hs.material.atomic_to_weight(list(comp.values()), list(comp.keys()))
            cout = {el: w/100 for el, w in zip(comp.keys(), weight_percent)}    # Divide by 100 to get percentage
            return cout

        self.get_comp = lambda x: _comp_weight_percent(x, _extract_info(params['composition']))
        #self.get_comp = _extract_info(params['composition'])
        if 'thickness' in params: self.get_thick = _extract_info(params['thickness'])
        if 'density' in params: self.get_dens = _extract_info(params['density'])

    def get_s(self, label=None):
        if label is None:
            return self.signals
        elif isinstance(label, str):
            for s in self.signals:
                if get_label(s) == label:
                    return s
        elif isinstance(label, list):
            s_out = []
            for s in self.signals:
                if get_label(s) in label:
                    s_out.append(s)
            return s_out
        else:
            raise Exception('The label (' + str(label) + ') has wrong datatype (' + str(type(label)) + '.')

    def _do_window_method(self):
        self.get_int = lambda x, lines: get_int_win(x, self.bw_pos.keys(), self.bw_pos, self.bw_width, self.iw_width)

    def print_int_mod(self, index):
        ints = self.get_int_mod(index)
        _print_ints(ints)

    def get_int_win(self, index):
        return get_int_win(self.signals[index], self.xray_lines, self._ext('bw_pos'), self._ext('bw_width'), self._ext('iw_width'))

    def print_int_win(self, index):
        ints = self.get_int_win(index)
        _print_ints(ints)

    def plot_window_results(self, index):
        return plot_window_results(self.signals[index], self.params['bw_pos'].keys(), self._ext('bw_pos'),
                                   self._ext('bw_width'), self._ext('iw_width'))

    def copy(self):
        return copy.deepcopy(self), copy.deepcopy(self.params)

    def get_signals(self, region):
        return sorted(list(filter(lambda x: get_region(x) == region, self.signals)), key=lambda x: get_ts(x))

    def save_report(self, size, *name):
        save_report(size, self.params['id']+'_', *name)

    def update(self, **kwargs):
        cp, cp_params = self.copy()
        cp_params.update(kwargs)
        cp.initialize(cp_params)
        return cp

    def _ext(self, key):
        if key in self.params:
            return self.params[key]
        else:
            return None

    def dump_obj(self, obj, label):
        with open(path_common + self.params['id'] + '_' + label + '.pkl', 'wb') as f:
            pickle.dump(obj, f)


class SpotAnalyzer(Analyzer):

    def __init__(self, params, calibration=None, create_tilt_series=True):
        self.initialize(params, calibration, create_tilt_series)

    def initialize(self, params, calibration=None, create_tilt_series=True):
        Analyzer.initialize(self, params)
        self.regions = self._ext('regions')
        if create_tilt_series:
            self.create_tilt_series()

    def _load_signals(self, params):
        self.folder_exports = self.folder + 'spectra\\'
        self.signals = load_files(self.folder_exports, '{s}{dd}', file_types=['ems', 'emsa'], signal_type='EDS_TEM',
                                  includeonly=self._ext('includeonly'))
        for s in self.signals:
            s.metadata.General.set_item('label', s.metadata.General.title)
            s.metadata.General.set_item('region', get_region(s))
            if 'region_names' in params:
                s.metadata.General.set_item('region_name', _extract_info(params['region_names'])(s))

    def create_tilt_series(self, regions=None, xray_lines=None):
        if regions is None: regions = self.regions
        if xray_lines is None: xray_lines = self.xray_lines
        self.tilt_series = create_tilt_series(self.signals, regions, xray_lines, self.get_int, self.get_comp,
                                         self.get_thick, self.get_dens, self.params['exclude'])
        for t in self.tilt_series:
            t.update({
                'zone_tilt': self._ext('zone_tilt'),
                'zone_axes': self._ext('zone_axes')
            })
        return self.tilt_series


    def dump_tilt_series(self):
        with open(path_common + self.params['id'] + '.pkl', 'wb') as f:
            pickle.dump(self.tilt_series, f)

    def plot_tilt_series(self, regions=None, xray_lines=None, avg=False,
                         annotate=False, ann_id=0, ann_off=-0.1, lines=True, zone_axes=True, zone_circles=False,
                         plot_exp_fit=False,
                         xmin=None, xmax=None, ymin=None, ymax=None, ref=None, legend=True, color='multi', ax_deg=True):
        return plot_tilt_series(self.tilt_series, regions, xray_lines, avg,
                                annotate, ann_id, ann_off, lines, zone_axes, zone_circles, plot_exp_fit,
                                xmin, xmax, ymin, ymax, ref, legend=legend, color=color, ax_deg=ax_deg)

    def get_zeta_factors(self, indices=None, regions=None, xray_lines=None, avg=None):
        return get_zeta_factors(self.tilt_series, indices, regions, xray_lines, avg)

    def plot_zeta_factors(self, indices=None, regions=None, xray_lines=None, avg=None):
        zf = get_zeta_factors(self.tilt_series, indices, regions, xray_lines, avg)
        return plot_zeta_factors(zf)

    def plot_signals(self, region):
        sgns = self.get_signals(region)
        hs.stack(sgns).plot(xray_lines='from_elements')
        for i, s in enumerate(sgns):
            print(str(i) + ': ' + s.metadata.General.title)

    def plot_effective_thickness(self, signals=True, thickness=100e-9, zone_tilt=True, zone_axes=True, beam_width=0.5):
        if zone_tilt is True:
            zone_tilt = self.params['zone_tilt']
        if zone_axes is True:
            zone_axes = self.params['zone_axes']
        if signals is True:
            signals = self.signals
        plot_effective_thickness(signals, thickness, zone_tilt, zone_axes, beam_width)

    def fit_shadowing(self, p0=[-20, 5, 2500], plot=False):
        return fit_shadowing(self.tilt_series, p0, plot)

    def quantify(self, label, method, factors):

        for s in self.signals:
            s.metadata.Sample.xray_lines = self.xray_lines  # Fix bug
            if 'beam_current' not in s.metadata.Acquisition_instrument.TEM:
                s.set_microscope_parameters(beam_current=1)    # Must be set, even though only necessary for thickness maps

        s = self.get_s(label)
        ints = [self.get_int[label][xray_line][0] for xray_line in self.xray_lines]
        if method == 'CL':
            qs = s.quantification(ints, 'CL', factors, composition_units='atomic')
        elif method == 'zeta':
            qs = s.quantification(ints, 'CL', factors, composition_units='atomic')
        return {q.metadata.Sample.xray_lines[0]: q.data[0] for q in qs}

    def plot_composition(self, kfacs, zfacs, signals=None):

        comps_cl = []
        comps_zeta = []
        labels = []

        for i, s in enumerate(self.get_s(signals)):
            label = get_label(s)
            labels.append(label)

            xray_line = self.xray_lines[0]

            comp_cl = self.quantify(label, 'CL', kfacs)
            comps_cl.append(comp_cl[xray_line])

            comp_zeta = self.quantify(label, 'zeta', zfacs)
            comps_zeta.append(comp_zeta[xray_line])

        fig, ax = plt.subplots()
        ax.plot(labels, comps_cl, 'o', label='Cliff-Lorimer')
        ax.plot(labels, comps_zeta, 'o', label='Zeta-factor method')
        plt.axhline(y=50, linestyle='--', color='lightgrey')
        plt.xticks(rotation=90)
        ax.set_ylabel(xray_line.split('_')[0] + ' composition (at%)')
        ax.set_xlabel('Spectrum ID')
        ax.legend()
        return fig, ax

    def plot_intensity_ratios(self, regions, line1, line2, annotate=False, legend=None):

        fig, ax = plt.subplots()

        ratios = {r: [] for r in regions}
        tilts = {r: [] for r in regions}
        labels = {r: [] for r in regions}

        for s in sorted(self.signals, key=lambda x: x.metadata.Acquisition_instrument.TEM.Stage.tilt_alpha):
            label = s.metadata.General.label
            region = get_region(label)
            if region not in regions:
                continue
            its = self.get_int[label]
            i1 = its[line1]
            i2 = its[line2]
            if i1 is None or i2 is None:
                continue
            ratio = i1[0].data[0] / i2[0].data[0]
            ratios[region].append(ratio)

            tilt = s.metadata.Acquisition_instrument.TEM.Stage.tilt_alpha
            tilts[region].append(tilt)

            labels[region].append(label)

        for i, region in enumerate(regions):
            tilt, ratio, label = tilts[region], ratios[region], labels[region]
            ax.plot(tilt, ratio, 'o', label=legend[i] if legend else None)
            if annotate:
                for l, t, r in zip(label, tilt, ratio):
                    ax.annotate(l, (t, r))

        if legend is not None:
            ax.legend()

        ax.set_xlabel('Stage tilt (degrees)')
        ax.set_ylabel('Ratio of I_' + _style_line(line1) + ' to I_' + _style_line(line2))


cache_pca = '.\\cache_pca\\'
cache_win_path = '.\\cache_win\\'


class MapAnalyzer(Analyzer):

    def initialize(self, params, calibration=None):
        Analyzer.initialize(self, params, calibration)

    def _load_signals(self, params):
        self.folder_exports = self.folder + 'maps\\'

        print('Loading signals...', end='')
        self.signals = load_files(self.folder_exports, '{s}{dd}', file_types=['hdf5', 'hspy'], signal_type='EDS_TEM',
                                  sort_by_ts=False, includeonly=self._ext('includeonly'))
        print('Done')
        for i, s in enumerate(self.signals):
            s.metadata.set_item('General.label', s.metadata.General.title)
            print('Changing datatype...', end='')
            s.change_dtype('float64')
            print('Done')
            if 'description' in params:
                s.metadata.set_item('General.description', _extract_info(params['description'])(s))

            if 'pre_processing' in params:
                print('Applying pre-processing step...', end='')
                s_proc = params['pre_processing'](s)
                self.signals[i] = s_proc
                print('Done')

        if params['pca'] is not None:
            comps = params['pca']
            s_out = []
            for s in self.signals:

                dim = len(s.axes_manager.shape)
                label = s.metadata.General.label + '_' + str(dim) + 'd_' + str(comps) + 'c.npz'
                if os.path.isdir(cache_pca) and label in os.listdir(cache_pca):
                    print('Retrieving PCA results for signal ' + s.metadata.General.label + ' from cache.')
                    s.learning_results.load(cache_pca + label)
                else:
                    print('PCA processing signal ' + s.metadata.General.label)
                    s.decomposition()
                    if not os.path.exists(cache_pca):
                        os.makedirs(cache_pca)
                    s.learning_results.save(cache_pca + label, overwrite=True)
                sd = s.get_decomposition_model(comps)
                sd.metadata.General.title = s.metadata.General.title
                s_out.append(sd)
            self.signals = s_out

    def _do_window_method(self):
        if not os.path.isdir(cache_win_path) or not self.params['intensity_cache']:
            self._create_win_int_cache(create_cache=self.params['intensity_cache'])
        else:
            self._get_win_int_from_cache()

    def _create_win_int_cache(self, create_cache=True):
        intensities = {}
        for s in self.signals:
            print('Calculating intensity of signal ' + s.metadata.General.label)
            intensity = get_int_win(s, list(self.bw_pos.keys()), self.bw_pos, self.bw_width, self.iw_width)
            intensities[get_label(s)] = {i.metadata.Sample.xray_lines[0]: [i] for i in intensity}

        if create_cache:
            if not os.path.exists(cache_win_path):
                os.makedirs(cache_win_path)
            for label, item in intensities.items():
                for line, intensity in item.items():
                    intensity[0].learning_results = LearningResults()    # In order to avoid HyperSpy bug
                    intensity[0].save(cache_win_path + label + '_' + line + '.hspy', overwrite=True)
        self.get_int = intensities

    def _get_win_int_from_cache(self):
        inp = {}
        for path in os.listdir(cache_win_path):
            s = hs.load(cache_win_path + path)
            print('Retrieving intensity of signal ' + s.metadata.General.label + ' from cache.')
            label = s.metadata.General.label
            xray_line = s.metadata.Sample.xray_lines[0]

            if label not in inp:
                inp[label] = {}
            inp[label][xray_line] = [s]
        self.get_int = inp

    def plot_zetas(self, labels=None, xray_lines=None, proc='auto', option='absolute', legend=True):
        return self._plot_something(labels, xray_lines, proc, option, legend)

    def plot_intensities(self, labels=None, xray_lines=None, proc='auto', legend=True):
        return self._plot_something(labels, xray_lines, proc, 'intensities', legend)

    def _plot_something(self, labels, xray_lines, proc, what, legend=True):
        """

        :param xray_lines:
        :param proc: Useful to be one of `lambda x: x.mean('y')` or `lambda x: x.inav[:, 0])`.
        """
        fig, ax = plt.subplots()
        if what == 'internal' or what == 'absolute':
            ys = self.get_zetas(labels, xray_lines, option=what, proc=proc)
        elif what == 'intensities':
            ys = self.get_intensities(labels, xray_lines, proc=proc)
        else:
            raise Exception('Unknown plot option in _plot_something: ' + what)
        for label, item in ys.items():
            s = self.get_s(label)
            for line, quantity in item.items():
                quantity = quantity[0]
                x = s.axes_manager.navigation_axes[0].axis
                unit = s.axes_manager.navigation_axes[0].units
                if legend == 'only_line': lab = _style_line(line)
                elif legend is True: lab = _style_line(line) + ', ' + s.metadata.General.description
                else: lab = None
                ax.plot(x, quantity, label=lab)
        if legend: ax.legend()
        ax.set_xlabel('Position (' + unit + ')')
        if what == 'absolute': ylabel = r'$\zeta$-factor (kg electron m$^{-2}$ photon$^{-1}$)'
        elif what == 'internal': ylabel = r'Internal $\zeta$-factor (photon$^{-1}$)'
        elif what == 'intensities': ylabel = 'Intensity (X-rays)'
        ax.set_ylabel(ylabel)
        ax.set_xlim(ax.get_xaxis().get_data_interval())
        return fig, ax

    def get_intensities(self, labels=None, xray_lines=None, proc='auto'):
        int_out = {}
        for label, item in self.get_int.items():
            if labels is None or label in labels:
                for line, signal in item.items():
                    if xray_lines is None or line in xray_lines:
                        if label not in int_out:
                            int_out[label] = {}
                        int_out[label][line] = [do_process(proc, signal[0])]
        if len(int_out) == 0:
            raise Exception('No elements found with specification: labels=' + str(labels), ', xray_lines=' + str(xray_lines))
        return int_out

    def get_intensity(self, label, xray_line, proc='auto'):
        return (list(list(self.get_intensities([label], [xray_line], proc).values())[0].values()))[0][0]

    def get_zetas(self, labels=None, xray_lines=None, option='absolute', proc='auto'):
        if xray_lines is None:
            xray_lines = self.xray_lines
        ints = self.get_intensities(labels, xray_lines, proc)
        zetas = {}
        for label, item in ints.items():
            s = self.get_s(label)
            if proc: s = do_process(proc, s)
            zetas[label] = {}
            for line, intensity in item.items():
                if option == 'internal': zeta = determine_relative_zeta_factor(s, [line], ints, self.get_comp)
                elif option == 'absolute': zeta = determine_zeta_factor(s, [line], ints, self.get_comp, self.get_thick, self.get_dens)
                else: raise Exception('Option ' + option + ' is not defined')
                zetas[label][line] = zeta
        return zetas

    def get_zeta(self, label, xray_line, proc):
        return (list(list(self.get_zetas([label], [xray_line], proc).values())[0].values()))[0][0]

    def _get_roi(self, s, x1, x2):
        xaxis = s.axes_manager[0].axis
        xrange = xaxis[-1] - xaxis[0]
        r = hs.roi.SpanROI(xaxis[0] + x1 * xrange, xaxis[0] + x2 * xrange)
        return r, r(s)

    def get_composition(self, label, method, factors, unknown=None, proc='auto'):

        xray_lines = list(factors.keys())
        facs = list(factors.values())
        if isinstance(facs[0], dict):
            facs = [fac['zeta'] for fac in facs]

        comps = {}
        if method == 'zeta_internal' or method == 'zeta_internal_normalized':
            for xray_line, fac in zip(xray_lines, facs):
                comps[xray_line] = np.float64(fac) * self.get_intensity(label, xray_line, proc)

            if unknown is not None:
                x = 1
                for xray_line, c in comps.items():
                    x = np.subtract(x, c)
                comps[unknown] = x

            if method == 'zeta_internal_normalized':
                comps_sum = np.zeros_like(list(comps.values())[0])
                for c in comps.values():
                    comps_sum += c
                for c in comps.values():
                    c /= comps_sum

            for xray_line, c in comps.items():
                comps[xray_line] = c*100    # # To get in percentage


        elif method == 'CL' or method == 'zeta':
            intensities = [self.get_int[label][xray_line][0] for xray_line in xray_lines]
            s = self.get_s(label)
            s.metadata.Sample.xray_lines = xray_lines
            if method == 'CL':    # k-factors specified by AZtec software are with respect to weight percentages
                qs = s.quantification(intensities, method='CL', factors=facs, composition_units='atomic')
            elif method == 'zeta':
                qs, tm = s.quantification(intensities, method='zeta', factors=facs, composition_units='atomic')
            comps = {q.metadata.Sample.xray_lines[0]: q for q in qs}

        else:
            raise Exception('Unknown quantification method: ' + method)

        return comps

    def quantify(self, label, method, factors):
        xray_lines = self.params['xray_lines']
        intensities = [self.get_int[label][xray_line][0] for xray_line in xray_lines]
        s = self.get_s(label)
        s.metadata.set_item('Sample.xray_lines', xray_lines)

        if method == 'CL':  # k-factors specified by AZtec software are with respect to weight percentages
            qs = s.quantification(intensities, method='CL', factors=factors, composition_units='atomic')
            return qs
        elif method == 'zeta':
            qs, tm = s.quantification(intensities, method='zeta', factors=factors, composition_units='atomic')
            return qs, tm

    def plot_composition(self, labels, method, factors, unknown=None, proc='auto', legend=True, broken_ylims=None, xray_lines=None):
        comp_list = []
        for label in labels:
            comp_list.append(self.get_composition(label, method, factors, unknown, proc))

        if broken_ylims is None:
            fig, ax = plt.subplots()
        else:
            from brokenaxes import brokenaxes
            fig = plt.figure()
            ax = brokenaxes(ylims=broken_ylims, d=None)
            ax2 = plt.gca() # hack
            ax2.spines['top'].set_visible(True)
            ax2.spines['right'].set_visible(True)
            ax.d = 0.015

        for i, comps in enumerate(comp_list):
            for xray_line, c in comps.items():
                if xray_lines is not None and xray_line not in xray_lines:
                    continue
                x = c.axes_manager.navigation_axes[0].axis
                unit = c.axes_manager.navigation_axes[0].units
                if legend == 'description': L = xray_line.split('_')[0] + ', ' + \
                                                self.get_s(labels[i]).metadata.General.description
                elif legend == 'element' or legend is True: L = xray_line.split('_')[0]
                else: L = None
                ax.plot(x, c, label=L)    # In percentage
        if legend: ax.legend()
        ax.set_xlabel('Position (' + unit + ')')
        ax.set_ylabel('Composition (at%)')

        if broken_ylims: ax.set_xlim(ax.get_xaxis()[0].get_data_interval())
        else: ax.set_xlim(ax.get_xaxis().get_data_interval())
        #ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda y, _: '{:.0%}'.format(y)))
        return fig, ax

    def determine_zeta_factor_from_internal(self, label, internal_zeta):
        if 'zeta' in internal_zeta:
            iz, xray_line = internal_zeta['zeta'], internal_zeta['xray_line']
        else:
            iz, xray_line = list(internal_zeta.keys())[0], list(internal_zeta.values())[0]

        zf = determine_zeta_factor_from_internal(self.get_s(label), [xray_line], lambda x: iz, self.get_thick,
                                                 self.get_dens)[0]
        return {'xray_line': xray_line, 'label': label, 'zeta': zf, 'std': None}

    def determine_zeta_factor(self, label, xray_line, rois=None, option='absolute', plot=False, plot_lines=None,
                              plot_label='only_line', comp=None):
        if option == 'internal':
            zf = determine_relative_zeta_factor(self.get_s(label), [xray_line], self.get_int,
                                                self.get_comp if comp is None else lambda x: comp)[0]
        elif option == 'absolute':
            zf = determine_zeta_factor(self.get_s(label), [xray_line], self.get_int,
                                       self.get_comp if comp is None else lambda x: comp, self.get_thick, self.get_dens)[0]
        else:
            raise Exception('Unknown option: ' + option)

        if rois is None: rois = [(0,1)]
        data = []
        roi_out = []
        for r in rois:
            roi, s = self._get_roi(zf, r[0], r[1])
            roi_out.append(roi)
            data.extend(list(np.reshape(s.data, -1)))

        zeta = np.mean(data)
        zeta_std = np.std(data)
        zdict = {'xray_line': xray_line, 'label': label, 'zeta': zeta, 'std': zeta_std}

        if plot:
            fig, ax = self.plot_intensities(labels=[label],
                                            xray_lines=plot_lines if plot_lines is not None else [xray_line],
                                            legend=plot_label)
            y_min, y_max = plt.ylim()
            for r in roi_out:
                rect = mpl.patches.Rectangle((r.left, y_min), width=r.right - r.left, height=y_max - y_min, alpha=0.3)
                ax.add_patch(rect)

        return zdict

    def plot_histogram(self, label, ks, zetas, xlim=(40, 60)):
        elem_id = 0
        q_zeta, mt = self.quantify(label, 'zeta', zetas)[elem_id]
        q_cl = self.quantify(label, 'CL', ks)[elem_id]
        fig, ax = plt.subplots()
        for q, mylabel in zip([q_cl, q_zeta], ['Cliff-Lorimer', 'Zeta-factors']):
            histo = q.get_histogram()
            x = histo.axes_manager['value'].axis
            ax.plot(x, histo.data, label=mylabel)
        ax.set_xlim(xlim)
        el = q.metadata.Sample.elements[0]
        ax.set_xlabel(el + ' composition (at%)')
        ax.set_ylabel('Frequency')
        ax.legend()
        plt.axvline(x=50, color='grey', linestyle='--')


def weight_to_atomic(weight_dict):
    atomic_percentages = hs.material.weight_to_atomic(list(weight_dict.values()))
    return {line: atp for line, atp in zip(list(weight_dict.keys()), atomic_percentages)}


def do_process(proc, signal):
    if proc=='auto':
        if len(signal.axes_manager.navigation_axes) == 2:
            proc = lambda x: x.sum('y')
        else:
            proc = None

    if proc is None:
        return signal
    else:
        return [proc(signal[0])]


def _print_ints(ints):
    for i in ints:
        print('{:1s}: {:2.2f}, \tratio: {:3.3f}'.format(i.metadata.Sample.xray_lines[0], i.data[0], i.data[0]/ints[0].data[0]))


def plot_effective_thickness(signals, thickness, zone_tilt, zone_axes, beam_width):
    fig, ax = plt.subplots()

    tilts = [i.metadata.Acquisition_instrument.TEM.Stage.tilt_alpha for i in signals]
    mn = min(tilts)
    mx = max(tilts)
    x = np.linspace(mn, mx, 1000)

    ax.plot(x, [get_thickness(i, thickness, zone_tilt) * 1e9 for i in x], label='Thin-film corrected', alpha=0.7)
    ax.plot(x, [get_thickness(i, thickness, zone_tilt, hexagonal=True) * 1e9 for i in x], alpha=0.7, label='Hexagonally corrected')
    ax.plot(x, [get_thickness(i, thickness, zone_tilt, hexagonal=True, beam_width=beam_width) * 1e9 for i in x], alpha=0.7, label='Beam width corrected')
    ax.legend()
    ax.set_ylabel('Sample thickness [nm]')
    ax.set_xlabel('Stage tilt')
    ax.xaxis.set_major_formatter(mpl.ticker.StrMethodFormatter(u'{x:.0f}°'))
    add_deviation(ax, min([get_thickness(i, thickness, zone_tilt)*1e9 for i in x]))
    plot_zone_axes(zone_axes)
    plt.tight_layout()


def get_thickness(s, width, zone_tilt, hexagonal=False, beam_width=0):
    if isinstance(s, hs.signals.BaseSignal): stage_tilt = s.metadata.Acquisition_instrument.TEM.Stage.tilt_alpha
    else: stage_tilt = s

    tilt = stage_tilt - zone_tilt
    if hexagonal:
        if tilt > 30:    # Degrees
            tilt = tilt - 60
    if beam_width == 0:
        fac = 1 / np.cos(np.deg2rad(tilt))
    else:    # Beam width correction (experimental)
        fac = (1 / np.cos(np.deg2rad(tilt)) - 1) * (1-beam_width) + 1
    t = get_uniform_thickness(width)
    return t*fac


def get_uniform_thickness(width):
    return np.sqrt(3)/2 * width


def load_ts(label, set_id=None):
    with open('../common/' + label + '.pkl', 'rb') as f:
        ts = pickle.load(f)
        if set_id is not None:
            for t in ts:
                t.update({'id': set_id})
        return ts


def func_shadowing(x, a, b, c):
    if isinstance(x, np.ndarray):
        return [_func_shadowing(i, a, b, c) for i in x]
    else:
        return _func_shadowing(x, a, b, c)


def _func_shadowing(x, a, b, c):
    if x >= b:
        return c
    elif x <= a:
        return np.infty
    else:
        return 1/((1/c)*(x-a)/(b-a))


def fit_shadowing(ts, regions=None, xray_lines=None, avg=False, p0=[-20, 5, 2500], plot=False):

    if regions is not None or xray_lines is not None or avg not in [False, None]:
        ts = ts_filter(ts, regions, xray_lines, avg)

        if len(ts) > 1:
            raise Exception('Filtering options cause more than one tilt series to be fitted.')
        elif len(ts) < 1:
            raise Exception('Filtering options cause no tilt series to be fitted.')
        ts = ts[0]

    res, cov = scipy.optimize.curve_fit(func_shadowing, ts['tilts'], ts['zetas'], p0=p0)

    if plot:
        plot_tilt_series([ts], legend=None, plot_exp_fit=True)

    return res, cov


def plot_shadowing(res, ts, legend=None, ax=None):
    xs = ts['tilts']
    mn, mx = xs[0], xs[-1]
    rg = mx - mn
    xfac = 0.02
    x = np.arange(mn - xfac * rg, mx + xfac * rg, 0.01)
    y = [func_shadowing(i, res[0], res[1], res[2]) for i in x]
    if ax is None:
        fix, ax = plt.subplots()
    ax.plot(x, y, label=legend)
    if legend:
        ax.legend()
