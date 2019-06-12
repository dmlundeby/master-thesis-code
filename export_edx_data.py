##################################################
# Script for exporting EDX maps from AZtec
##################################################

print('Loading Python libraries')

import os
import sys
import hyperspy.api as hs
import numpy as np
import json

print()

path = '.\\'

for arg in sys.argv[1:]:
    if not arg.startswith('-'):
        path = arg.rstrip('\\') + '\\'
        break


transpose = False
if '--transpose' in sys.argv:
    transpose = True

flipx = False
if '--flipx' in sys.argv:
    flipx = True

flipy = False
if '--flipy' in sys.argv:
    flipy = True

def main():
    labels = []
    contents = os.listdir(path)
    for filename in contents:
        fn_spl = filename.split('.')
        if len(fn_spl) == 1:
            continue
        label, filetype = ''.join(fn_spl[:-1]), fn_spl[-1]
        if filetype == 'raw':
            labels.append(label)

    for label in labels:
        print('Processing', label)
        pathlabel = path + '\\' + label
        convert_aztec_to_hdf5(pathlabel, transpose, flipx, flipy)

    if len(labels) == 0:
        print('Couldn\'t find any files to process')


def convert_aztec_to_hdf5(pathlabel, transpose=False, flipx=False, flipy=False):
    """
    Convert AZtec EDX map referenced by `pathlabel` to an `.hdf5` file, and save this file to the same folder.

    The function assumes the following files exist:

    * `.raw`: The raw EDX data.
    * `.rpl`: Information about the structure of the raw data.
    * `.emsa`: Sum spectrum, with calibration data and metadata.
    * `.txt`: Additional metadata.

    :param pathlabel: Full path of the files containing EDX information, including label but without extension. Example: `C:/Users/john/Documents/EDX/A1`. In this case, `A1` is the label and the following files are assumed to exist:  `A1.raw`, `A1.rpl`, `A1.emsa`, `A1.txt`.
    :param transpose: Whether to interchange the `x` and `y` axes.
    :param flipx: Whether to make the `x` axis reversed.
    :param flipy: Whether to make the `y` axis reversed.
    :return: The `.hdf5`signal which has been created.
    """

    fn_rpl = pathlabel + '.rpl'
    fn_emsa = pathlabel + '.emsa'
    fn_txt = pathlabel + '.txt'


    # Load spectrum
    # The .rpl file uses data from the .raw file.
    # For some reason, the signal needs to be transposed.
    s = hs.load(fn_rpl, signal_type='EDS_TEM').T

    # Load sum spectrum
    # The sum spectrum contains some metadata which is not contained in the .rpl file
    s_sum = hs.load(fn_emsa, signal_type='EDS_TEM')
    s.get_calibration_from(s_sum)

    # The EDX detector real time and live time needs to divided by the number of pixels
    pixels=1
    for axis in s.axes_manager.navigation_axes:
        pixels *= axis.size
    s.metadata.Acquisition_instrument.TEM.Detector.EDS.real_time /= pixels
    s.metadata.Acquisition_instrument.TEM.Detector.EDS.live_time /= pixels

    # Some properties are not transfered by s.get_calibration_from
    def transfer(item):
        s.metadata.set_item(item, s_sum.metadata.get_item(item))
    transfer('General.date')
    transfer('General.time')

    # Load metadata
    meta_np = np.genfromtxt(fn_txt, encoding='utf-8-sig', delimiter=':\t', dtype='str')
    meta = dict(meta_np)

    # Set properties from metadata
    ax_ids = ['width', 'height']
    ax_names = ['x', 'y']
    sizes_id = ['Image Width', 'Image Height']
    pxs_id = ['Resolution (Width)', 'Resolution (Height)']

    for ax_id, ax_name, size_id, px_id in zip(ax_ids, ax_names, sizes_id, pxs_id):
        ax = s.axes_manager[ax_id]
        size, size_unit = split_property(meta[size_id])
        px, px_unit = split_property(meta[px_id])
        ax.scale = float(size.replace(',', '.')) / float(px)
        ax.units = size_unit
        ax.name = ax_name

    if transpose:
        s = s.transpose(signal_axes=[2], navigation_axes=[1, 0])
        s.axes_manager[0].name = ax_names[0]
        s.axes_manager[1].name = ax_names[1]

    if flipx:
        s.data = np.flip(s.data, axis=1)

    if flipy:
        s.data = np.flip(s.data, axis=0)

    s.save(pathlabel, extension='hdf5', overwrite=True)

    return s


def split_property(text):
    lim = None
    for i in range(len(text) - 1, -1, -1):
        if text[i].isnumeric():
            lim = i + 1
            break
    return text[:lim], text[lim:]


def get_index(string):
    s = ".".join(string.split('.')[:-1]).split('-')[0]
    start = -1
    for i in range(len(s)-1, -1, -1):
        if not s[i].isdigit():
            start = i
            break
    return s[start+1:]


def load_json(path):
    with open(path) as file:
        return json.load(file)

main()

input('Press Enter to exit')
