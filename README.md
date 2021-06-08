# ADCIRCPy

## Python library for automating ADCIRC model runs.

[![tests](https://github.com/JaimeCalzadaNOAA/adcircpy/workflows/tests/badge.svg)](https://github.com/JaimeCalzadaNOAA/adcircpy/actions?query=workflow%3Atests)
[![build](https://github.com/JaimeCalzadaNOAA/adcircpy/workflows/build/badge.svg)](https://github.com/JaimeCalzadaNOAA/adcircpy/actions?query=workflow%3Abuild)
[![codecov](https://codecov.io/gh/JaimeCalzadaNOAA/adcircpy/branch/main/graph/badge.svg?token=BQWB1QKJ3Q)](https://codecov.io/gh/JaimeCalzadaNOAA/adcircpy)
[![version](https://img.shields.io/pypi/v/adcircpy)](https://pypi.org/project/adcircpy)
[![license](https://img.shields.io/github/license/JaimeCalzadaNOAA/adcircpy)](https://opensource.org/licenses/gpl-license)
[![style](https://sourceforge.net/p/oitnb/code/ci/default/tree/_doc/_static/oitnb.svg?format=raw)](https://sourceforge.net/p/oitnb/code)

### Installation notes:

Please use a virtual environment with Python>=3.6. You may use conda or the OS's Python to provide a virtual environment for
the application.

You may install the application though pip. This will install the latest tagged version.
<!--pytest-codeblocks:skip-->

```bash
pip install adcircpy
```

Alternatively, you many manually install the repo by cloning it and then running
<!--pytest-codeblocks:skip-->

```bash
pip install .
```

### Command Line:

This program exposes a few commands available from the command line interface. You may pass the `-h` flag to any of this
commands to explore their functionality.

* `tide_gen`
* `plot_mesh`
* `tidal_run`
* `best_track_run`
* `best_track_file`
* `plot_maxele`
* `plot_fort61`
* `fort63`

#### Command line examples:

##### Generate tidal constituent template from command line

You can quickly create a tidal component table for your your mesh by executing the `tide_gen` command and by passing a mesh, a
start date and number of run days as arguments. This functions sources data from
the [HAMTIDE](https://icdc.cen.uni-hamburg.de/en/hamtide.html) database by default.
<!--pytest-codeblocks:skip-->

```bash
tide_gen \
    /path/to/your/fort.14 \
    2021-02-26T00:00:00 \
    15 \
    --mesh-crs='epsg:4326'
```

##### Hurricane Sandy (AL182012)

To create the ADCIRC input files includes both tides and storm data for Hurricane Sandy:
<!--pytest-codeblocks:skip-->

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

Note that the --crs flag is required due to the fort.14 not containing Coordinate Reference System information which is
required for correct operation. [EPSG:4326](https://spatialreference.org/ref/epsg/wgs-84/) means that the mesh is in WGS84 (
lat/lon). Note that the backlash represents "continue on next line" for the shell. You may write the command above on a single
line after excluding the backslashes.

##### Quick plots

These are two examples for doing quick plots with the package. These are given here as illustrative examples only. There is
support for more file types than this examples, but the program does not yet support every output input/output file type. As a
user, you are encouraged to explore what's available and suggest and contribute your improvements.
<!--pytest-codeblocks:skip-->

```bash
plot_fort61 /path/to/fort.61.nc MSL --show --coops-only
```

<!--pytest-codeblocks:skip-->

```bash
plot_mesh /path/to/fort.14 --show-elements
```

### Python API:

See the [examples](examples) directory for usage examples.

#### `example_1.py`

The following code builds a simple ADCIRC run configuration by doing the following:

1. reads a `fort.14` mesh file (specifically a test mesh for Shinnecock Inlet)
2. adds tidal forcing to the mesh
3. creates an `AdcircRun` driver object with the mesh, including start and end dates
4. overrides default model options in the resulting `fort.15`
5. runs ADCIRC if present, otherwise writes configuration to disk

```python
from datetime import datetime, timedelta
import logging
from pathlib import Path
import shutil

from adcircpy import AdcircMesh, AdcircRun, Tides
from adcircpy.utilities import download_mesh

MESH_URL = 'https://www.dropbox.com/s/1wk91r67cacf132/NetCDF_shinnecock_inlet.tar.bz2?dl=1'

DATA_DIRECTORY = Path(__file__).parent.absolute() / 'data'
INPUT_DIRECTORY = DATA_DIRECTORY / 'input' / 'NetCDF_Shinnecock_Inlet'
OUTPUT_DIRECTORY = DATA_DIRECTORY / 'output' / 'example_1'

download_mesh(
    url=MESH_URL, directory=INPUT_DIRECTORY,
)

# open mesh file
mesh = AdcircMesh.open(INPUT_DIRECTORY / 'fort.14', crs=4326)

# initialize tidal forcing and constituents
tidal_forcing = Tides()
tidal_forcing.use_constituent('M2')
tidal_forcing.use_constituent('N2')
tidal_forcing.use_constituent('S2')
tidal_forcing.use_constituent('K1')
tidal_forcing.use_constituent('O1')
mesh.add_forcing(tidal_forcing)

# set simulation dates
duration = timedelta(days=5)
start_date = datetime(2015, 12, 14)
end_date = start_date + duration

# instantiate driver object
driver = AdcircRun(mesh, start_date, end_date)

# request outputs
driver.set_elevation_surface_output(sampling_rate=timedelta(minutes=30))
driver.set_velocity_surface_output(sampling_rate=timedelta(minutes=30))

# override default options so the resulting `fort.15` matches the original Shinnecock test case options
driver.timestep = 6.0
driver.DRAMP = 2.0
driver.TOUTGE = 3.8
driver.TOUTGV = 3.8
driver.smagorinsky = False
driver.horizontal_mixing_coefficient = 5.0
driver.gwce_solution_scheme = 'semi-implicit-legacy'

if shutil.which('padcirc') is not None:
    driver.run(OUTPUT_DIRECTORY, overwrite=True)
elif shutil.which('adcirc') is not None:
    driver.run(OUTPUT_DIRECTORY, overwrite=True, nproc=1)
else:
    logging.warning(
        'ADCIRC binaries were not found in PATH. '
        'ADCIRC will not run. Writing files to disk...'
    )
    driver.write(OUTPUT_DIRECTORY, overwrite=True)
```

#### `example_2.py`

The following code is similar to `example_1.py`, above, except it adds a static Manning's N coefficient to the mesh.

```python
from datetime import datetime, timedelta
import logging
from pathlib import Path
import shutil

import numpy

from adcircpy import AdcircMesh, AdcircRun, Tides
from adcircpy.utilities import download_mesh

MESH_URL = 'https://www.dropbox.com/s/1wk91r67cacf132/NetCDF_shinnecock_inlet.tar.bz2?dl=1'

DATA_DIRECTORY = Path(__file__).parent.absolute() / 'data'
INPUT_DIRECTORY = DATA_DIRECTORY / 'input' / 'NetCDF_Shinnecock_Inlet'
OUTPUT_DIRECTORY = DATA_DIRECTORY / 'output' / 'example_2'

download_mesh(
    url=MESH_URL, directory=INPUT_DIRECTORY,
)

# open mesh file
mesh = AdcircMesh.open(INPUT_DIRECTORY / 'fort.14', crs=4326)

# generate tau0 factor
mesh.generate_tau0()

# also add Manning's N to the domain (constant for this example)
mesh.mannings_n_at_sea_floor = numpy.full(mesh.values.shape, 0.025)

# initialize tidal forcing and constituents
tidal_forcing = Tides()
tidal_forcing.use_constituent('M2')
tidal_forcing.use_constituent('N2')
tidal_forcing.use_constituent('S2')
tidal_forcing.use_constituent('K1')
tidal_forcing.use_constituent('O1')
mesh.add_forcing(tidal_forcing)

# set simulation dates
spinup_time = timedelta(days=2)
duration = timedelta(days=3)
start_date = datetime(2015, 12, 14) + spinup_time
end_date = start_date + duration

# instantiate driver object
driver = AdcircRun(mesh, start_date, end_date, spinup_time)

# request outputs
driver.set_elevation_surface_output(sampling_rate=timedelta(minutes=30))
driver.set_velocity_surface_output(sampling_rate=timedelta(minutes=30))

# override default options
driver.timestep = 4.0

if shutil.which('padcirc') is not None:
    driver.run(OUTPUT_DIRECTORY, overwrite=True)
elif shutil.which('adcirc') is not None:
    driver.run(OUTPUT_DIRECTORY, overwrite=True, nproc=1)
else:
    logging.warning(
        'ADCIRC binaries were not found in PATH. '
        'ADCIRC will not run. Writing files to disk...'
    )
    driver.write(OUTPUT_DIRECTORY, overwrite=True)
```

#### `example_3.py`

The following code is similar to `example_1.py`, above, except it adds HURDAT BestTrack wind forcing and also builds a Slurm
job script for submission to a job manager.

```python
from datetime import timedelta
from pathlib import Path

from adcircpy import AdcircMesh, AdcircRun, Tides
from adcircpy.forcing.winds import BestTrackForcing
from adcircpy.server import SlurmConfig
from tests import download_mesh

MESH_URL = 'https://www.dropbox.com/s/1wk91r67cacf132/NetCDF_shinnecock_inlet.tar.bz2?dl=1'

DATA_DIRECTORY = Path(__file__).parent.absolute() / 'data'
INPUT_DIRECTORY = DATA_DIRECTORY / 'input' / 'NetCDF_Shinnecock_Inlet'
OUTPUT_DIRECTORY = DATA_DIRECTORY / 'output' / 'example_3'

download_mesh(
    url=MESH_URL, directory=INPUT_DIRECTORY,
)

# open mesh file
mesh = AdcircMesh.open(INPUT_DIRECTORY / 'fort.14', crs=4326)

# initialize tidal forcing and constituents
tidal_forcing = Tides()
tidal_forcing.use_all()
mesh.add_forcing(tidal_forcing)

# initialize wind forcing
wind_forcing = BestTrackForcing('Sandy2012')
mesh.add_forcing(wind_forcing)

# initialize Slurm configuration
slurm = SlurmConfig(
    account='account',
    ntasks=1000,
    run_name='adcircpy/examples/example_3.py',
    partition='partition',
    walltime=timedelta(hours=8),
    mail_type='all',
    mail_user='example@email.gov',
    log_filename='example_3.log',
    modules=['intel/2020', 'impi/2020', 'netcdf/4.7.2-parallel'],
    path_prefix='$HOME/adcirc/build',
)

# instantiate driver object
driver = AdcircRun(mesh, spinup_time=timedelta(days=15), server_config=slurm)

# write driver state to disk
driver.write(OUTPUT_DIRECTORY, overwrite=True)
```

#### `example_4.py`

The following code is similar to `example_3.py`, above, except it uses ATMESH wind forcing and WW3DATA wave forcing.

```python
from datetime import datetime, timedelta
from pathlib import Path

from adcircpy import AdcircMesh, AdcircRun, Tides
from adcircpy.forcing.waves.ww3 import WaveWatch3DataForcing
from adcircpy.forcing.winds.atmesh import AtmosphericMeshForcing
from adcircpy.server import SlurmConfig
from tests import download_mesh

MESH_URL = 'https://www.dropbox.com/s/1wk91r67cacf132/NetCDF_shinnecock_inlet.tar.bz2?dl=1'

DATA_DIRECTORY = Path(__file__).parent.absolute() / 'data'
INPUT_DIRECTORY = DATA_DIRECTORY / 'input' / 'NetCDF_Shinnecock_Inlet'
OUTPUT_DIRECTORY = DATA_DIRECTORY / 'output' / 'example_4'

download_mesh(
    url=MESH_URL, directory=INPUT_DIRECTORY,
)

# open mesh file
mesh = AdcircMesh.open(INPUT_DIRECTORY / 'fort.14', crs=4326)

# initialize tidal forcing and constituents
tidal_forcing = Tides()
tidal_forcing.use_all()
mesh.add_forcing(tidal_forcing)

# initialize atmospheric mesh forcings (for NUOPC coupling)
wind_forcing = AtmosphericMeshForcing(
    filename='Wind_HWRF_SANDY_Nov2018_ExtendedSmoothT.nc', nws=17, interval_seconds=3600,
)
mesh.add_forcing(wind_forcing)

# initialize wave mesh forcings (for NUOPC coupling)
wave_forcing = WaveWatch3DataForcing(
    filename='ww3.HWRF.NOV2018.2012_sxy.nc', nrs=5, interval_seconds=3600,
)
mesh.add_forcing(wave_forcing)

# initialize Slurm configuration
slurm = SlurmConfig(
    account='account',
    ntasks=1000,
    run_name='adcircpy/examples/example_4.py',
    partition='partition',
    walltime=timedelta(hours=8),
    mail_type='all',
    mail_user='example@email.gov',
    log_filename='example_4.log',
    modules=['intel/2020', 'impi/2020', 'netcdf/4.7.2-parallel'],
    path_prefix='$HOME/adcirc/build',
)

# instantiate driver object
driver = AdcircRun(
    mesh=mesh,
    start_date=datetime.now(),
    end_date=timedelta(days=7),
    spinup_time=timedelta(days=5),
    server_config=slurm,
)

# write driver state to disk
driver.write(OUTPUT_DIRECTORY, overwrite=True)
```

### Contact

For questions comments and suggestions, please email me at [jreniel@gmail.com](mailto:jreniel@gmail.com).
