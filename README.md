# AdcircPy

### A Python interface for handling inputs and outputs for the ADCIRC hydrodynamic model.

### Basic Installation:

You can install this software through conda or through pip.</br>
The dependecies of Python packages installed through pip must be satisfied by the operating system (see Method 2 for details). On a normal user's computer, the recommended method of installation is through conda. Whether you are using the pip method or the conda method, it is generally good practice to always make use of a virtual environment for your project. If you are a Python newcomer, make sure to read upon Python environments and their usage as you dive more into Python development.</br>

If you are doing development on a Windows machine it is highly recommended to use the Windows Subsystem for Linux, and run this software from a Widows-Subsystem-for-Linux shell. MingW64, Cygwin, the git-shell-for-windows, or any other off-the-shelf POSIX layer on Windows should work.</br>

### Method 1: Conda

The advantage of using conda is that it provides the necessary dependencies as precompiled binaries, therefore no local compilations are required, and cuts the installation time substantially. This is the recommended method to use. </br>

Install Miniconda3 (recommended) or Anaconda3 on your system. Rememeber to make this installation on the POSIX layer when using Windows.</br>

Create a new conda environment for your project, and inside the environment run:

```cmd
conda update -n base -c defaults conda
conda install -c conda-forge gdal
pip install AdcircPy
```

### Method 2: Pip-only (requires compilation of system libraries).

Pip requires that some system libraries are precompiled and installed before running pip install command. _It is highly recommended that you use the system packge manager to satisfy these dependencies_, or, they may be compiled from source. The [Dockerfile](./Dockerfile) contains a full build on an Alpine distribution, which may be used for reference. Note that "dev" versions mean that the headers are required along with the compiled library. Also note that some of these libraries you may already have in your system.

#### Full system dependency list:

- g++
- gcc
- m4
- make
- autoconf
- zlib-dev
- libc-dev
- curl-dev
- linux-headers
- hdf5-dev
- netcdf-c-dev
- python3-dev
- sqlite-dev
- Proj4-dev
- GDAL-dev
- openblas-dev (for scipy)
- freetype-dev (for matplotlib)
- py3-pip

Once these dependencies have been met, you may simply do

```cmd
pip install AdcircPy
```

### Ugrading:

If you have installed AdcircPy and you would like to upgrade, simply do:

```cmd
pip install AdcircPy -U
```

You can always run this command to make sure you have the latest version since the latest in PyPi should correspond to the latest tagged version of this repository.

#### NOTE: If you are developer and would like to debug the package, you may install this package on developer mode by first cloning the package, and then doing:

```cmd
pip install -e .
```

### Usage

#### CLI interface

After installation a few new commands are available on the command line interface.

A few of the available commands are:

- PlotMesh
- GenerateTidalRun
- GenerateBestTrackFile
- PlotMaxele
- HighWaterMarkValidation
- PlotTidalStationsOutput

Since this package is in development, the list above may not be up to date, therefore you may check [setup.py](setup.py), where the entrypoints for the present release are listed.

#### As an API

This package is an API to handle input and output files for the [ADCIRC](http://adcirc.org) hydrodynamic model.
You can load an Adcirc mesh and plot it by doing:

```Python
from adcircpy.mesh import AdcircMesh
mesh = AdcircMesh.open('/path/to/fort.14', SpatialReference=4326)
# The most basic you can do with this software is create plots:
mesh.make_plot(show=True)
```

### Example of fort.15 generation:
```Python
#! /usr/bin/env python
"""
Isabel 2003 - DelawareBayOceanMesh2D
"""
import os
from datetime import timedelta
from adcircpy.mesh import AdcircMesh
from adcircpy.model import BestTrackForcing, TidalForcing, AdcircRun


def main():
    #  ----------- instantiante run
    adcirc_run = AdcircRun()
    adcirc_run.AdcircMesh = AdcircMesh.open(
        os.path.abspath('../fort.14'), 4326)

    #  ----------- add nodal attributes
    adcirc_run.AdcircMesh.import_nodal_attributes(
        os.path.abspath('../fort.13'))
    for attribute in adcirc_run.AdcircMesh.get_nodal_attribute_names():
        adcirc_run.AdcircMesh.set_nodal_attribute_state(attribute, True, True)

    #  ----------- add tidal forcing
    adcirc_run.TidalForcing = TidalForcing()
    adcirc_run.TidalForcing.use_major()

    #  ----------- add wind forcing
    adcirc_run.WindForcing = BestTrackForcing()
    adcirc_run.WindForcing.storm_id = 'AL182012'
    adcirc_run.WindForcing.remove_TS()
    # adcirc_run.WindForcing.remove_EX()

    #  ----------- request outputs
    adcirc_run.copy_fort15_stations(os.path.abspath('stations.txt'))
    adcirc_run.set_elevation_stations_output(timedelta(minutes=6.))
    adcirc_run.set_elevation_global_output(timedelta(minutes=60.))

    #  ----------- set additional run options
    adcirc_run.start_date = adcirc_run.WindForcing.start_date
    adcirc_run.end_date = adcirc_run.WindForcing.end_date
    adcirc_run.spinup_time = timedelta(days=15.)
    adcirc_run.gwce_solution_scheme = 'explicit'
    adcirc_run.DTDP = 2.

    #  ----------- dump run
    os.makedirs(os.path.abspath('./'), exist_ok=True)
    adcirc_run.dump(os.path.abspath('./'), overwrite=True)


if __name__ == '__main__':
    main()
```
