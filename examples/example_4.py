from datetime import datetime, timedelta
from pathlib import Path

from adcircpy import AdcircMesh, AdcircRun, Tides
from adcircpy.forcing.waves.ww3 import WaveWatch3DataForcing
from adcircpy.forcing.winds.atmesh import AtmosphericMeshForcing
from adcircpy.server import SlurmConfig
from adcircpy.utilities import download_mesh

DATA_DIRECTORY = Path(__file__).parent.absolute() / 'data'
INPUT_DIRECTORY = DATA_DIRECTORY / 'input'
OUTPUT_DIRECTORY = DATA_DIRECTORY / 'output' / 'example_4'

MESH_DIRECTORY = INPUT_DIRECTORY / 'shinnecock'

download_mesh(
    url='https://www.dropbox.com/s/1wk91r67cacf132/NetCDF_shinnecock_inlet.tar.bz2?dl=1',
    directory=MESH_DIRECTORY,
    known_hash='99d764541983bfee60d4176af48ed803d427dea61243fa22d3f4003ebcec98f4',
)

# open mesh file
mesh = AdcircMesh.open(MESH_DIRECTORY / 'fort.14', crs=4326)

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
