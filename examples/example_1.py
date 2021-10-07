from datetime import datetime, timedelta
from pathlib import Path
import shutil

from adcircpy import AdcircMesh, AdcircRun, Tides
from adcircpy.utilities import download_mesh, get_logger

LOGGER = get_logger(__name__)

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
    LOGGER.warning(
        'ADCIRC binaries were not found in PATH. '
        'ADCIRC will not run. Writing files to disk...'
    )
    driver.write(OUTPUT_DIRECTORY, overwrite=True)
