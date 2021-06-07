from datetime import datetime, timedelta
from pathlib import Path
import shutil
import tempfile

from adcircpy.driver import AdcircRun
from adcircpy.forcing import Tides
from adcircpy.mesh import AdcircMesh
from tests import download_mesh

DATA_DIRECTORY = Path(__file__).parent.absolute().resolve() / 'data'
INPUT_DIRECTORY = DATA_DIRECTORY / 'input' / 'NetCDF_Shinnecock_Inlet'

download_mesh(
    url='https://www.dropbox.com/s/1wk91r67cacf132/NetCDF_shinnecock_inlet.tar.bz2?dl=1',
    directory=INPUT_DIRECTORY,
)


def test_tidal_run():
    mesh = AdcircMesh.open(INPUT_DIRECTORY / 'fort.14', crs=4326)

    tidal_forcing = Tides()
    tidal_forcing.use_all()

    mesh.add_forcing(tidal_forcing)
    now = datetime.utcnow()
    driver = AdcircRun(
        mesh,
        start_date=now,
        end_date=now + timedelta(days=0.5),
        spinup_time=timedelta(days=0.5),
    )
    driver.timestep = 10.0

    if shutil.which('padcirc') is None:
        tmpdir = tempfile.TemporaryDirectory()
        driver.write(tmpdir.name)
    else:
        driver.run()
