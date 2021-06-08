from datetime import datetime, timedelta
from pathlib import Path

from adcircpy import AdcircMesh, AdcircRun, Tides
from adcircpy.forcing.winds import BestTrackForcing
from adcircpy.server import SlurmConfig
from adcircpy.utilities import download_mesh

DATA_DIRECTORY = Path(__file__).parent.absolute() / 'data'
INPUT_DIRECTORY = DATA_DIRECTORY / 'input'
OUTPUT_DIRECTORY = DATA_DIRECTORY / 'output' / 'example_3'

MESH_URL = 'https://www.dropbox.com/s/1wk91r67cacf132/NetCDF_shinnecock_inlet.tar.bz2?dl=1'
MESH_DIRECTORY = INPUT_DIRECTORY / 'shinnecock'

download_mesh(
    url=MESH_URL, directory=MESH_DIRECTORY,
)

# open mesh file
mesh = AdcircMesh.open(MESH_DIRECTORY / 'fort.14', crs=4326)

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

# set simulation dates
spinup_time = timedelta(days=15)
duration = timedelta(days=3)
start_date = datetime(2012, 10, 21, 18)
end_date = start_date + duration

# instantiate driver object
driver = AdcircRun(mesh, start_date, end_date, spinup_time, server_config=slurm)

# write driver state to disk
driver.write(OUTPUT_DIRECTORY, overwrite=True)
