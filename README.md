# Code written as part of my Master's thesis

Title: "Improving the accuracy of TEM-EDX quantification by implementing the &#950;-factor method".

The project has involved development of code for EDX data analysis in Python, as well and some scripts.

## Python code for analyzing EDX data

`uts.py`: Python module containing classes and functions for analyzing EDX data. The module is centered around two main classes: `SpotAnalyzer` (for analyzing single EDX spectra) and `MapAnalyzer` (for analyzing EDX maps).


## DigitalMicrograph scripts

`find_beam_current.s`: Script for measuring the TEM beam current using the CCD camera.

`save_metadata.s`: Script for saving metadata from the TEM (such as the stage tilt) to a `.json` file. This metadata is not saved otherwise.


## Code for exporting EDX data from AZtec software

`export_edx_data.py`: Python script for exporting EDX maps from the AZtec software to `.hdf5` format. The script requires `.raw` and `.rpl` files as well as some metatadata files as input:

`export_edx_data.au3`: Simple Windows utility (written in AutoIt) which acts as a wrapper to the `export_edx_data.py` script.

`update_metadata.py`: Script which updates metadata of EDX spectra, using `.json` files from the `save_metadata.s` script.


## Jupyter notebooks

`calculating_self_absorption_plots.ipynb`: Code for creating a plot to visualize self absorption, using absorption coefficiencts from HyperSpy.
