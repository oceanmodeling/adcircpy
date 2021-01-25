# ADCIRCPy
## Python library for automating ADCIRC model runs.
[![tests](https://github.com/JaimeCalzadaNOAA/adcircpy/workflows/tests/badge.svg)](https://github.com/JaimeCalzadaNOAA/adcircpy/actions?query=workflow%3Atests)
[![build](https://github.com/JaimeCalzadaNOAA/adcircpy/workflows/build/badge.svg)](https://github.com/JaimeCalzadaNOAA/adcircpy/actions?query=workflow%3Abuild)
[![coverage](tests/coverage.svg)](https://github.com/JaimeCalzadaNOAA/adcircpy/actions)
[![version](https://img.shields.io/pypi/v/adcircpy)](https://pypi.org/project/adcircpy)
[![license](https://img.shields.io/github/license/JaimeCalzadaNOAA/adcircpy)](https://opensource.org/licenses/gpl-license)

### Installation notes:

Please use a virtual environment with Python>=3.6. You may use conda or the OS's Python to provide a virtual environment for the application.

You may install the application though pip. This will install the latest tagged version.
```bash
pip install adcircpy
```

Alternatively, you many manually install the repo by cloning it and then running
```bash
pip install .
```

### Examples: 
See the [examples](examples) directory for usage examples.


### Command Line:
This program exposes a few commands available from the command line interface. You may pass the `-h` flag to any of this commands to explore their functionality. 
* `plot_mesh`
* `tidal_run`
* `best_track_run`
* `best_track_file`
* `plot_maxele`
* `plot_fort61` 
* `fort63`

#### Command line examples:
##### Hurricane Sandy (AL182012)
To create the ADCIRC input files includes both tides and storm data for Hurricane Sandy:
```bash
best_track_run \
    /path/to/your/fort.14 \
    Sandy2012 \
    --fort13=/path/to/your/fort.13 \
    --crs=EPSG:4326 \
    --output-directory=/path/where/you/want/the/files \
    --constituents=all \
    --spinup-days=15.0 \
    --elev=30. \
    --mete=30. \
    --velo=30. \
    --skip-run
```
Note that the --crs flag is required due to the fort.14 not containing Coordinate Reference System information which is required for correct operation. [EPSG:4326](https://spatialreference.org/ref/epsg/wgs-84/) means that the mesh is in WGS84 (lat/lon).
Note that the backlash represents "continue on next line" for the shell. You may write the command above on a single line after excluding the backslashes.

##### Quick plots
These are two examples for doing quick plots with the package. These are given here as illustrative examples only. There is support for more file types than this examples, but the program does not yet support every output input/output file type.
As a user, you are encouraged to explore what's available and suggest and contribute your improvements.
```bash
plot_fort61 /path/to/fort.61.nc MSL --show --coops-only
```
```bash
plot_mesh /path/to/fort.14 --show-elements
```

### Contact
For questions comments and suggestions, please email me at jreniel@gmail.com
