from datetime import datetime, timedelta
from pathlib import Path
import shutil
import warnings

import numpy

from adcircpy import AdcircMesh, AdcircRun, Tides
from adcircpy.utilities import download_mesh

DATA_DIRECTORY = Path(__file__).parent.absolute() / 'data'
INPUT_DIRECTORY = DATA_DIRECTORY / 'input'
OUTPUT_DIRECTORY = DATA_DIRECTORY / 'output' / 'example_2'

MESH_DIRECTORY = INPUT_DIRECTORY / 'shinnecock'

download_mesh(
    url='https://www.dropbox.com/s/1wk91r67cacf132/NetCDF_shinnecock_inlet.tar.bz2?dl=1',
    directory=MESH_DIRECTORY,
    known_hash='99d764541983bfee60d4176af48ed803d427dea61243fa22d3f4003ebcec98f4',
)

# open mesh file
mesh = AdcircMesh.open(MESH_DIRECTORY / 'fort.14', crs=4326)

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
    warnings.warn(
        'ADCIRC binaries were not found in PATH. '
        'ADCIRC will not run. Writing files to disk...'
    )
    driver.write(OUTPUT_DIRECTORY, overwrite=True)
