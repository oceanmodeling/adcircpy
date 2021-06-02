from datetime import datetime, timedelta
from pathlib import Path
import tempfile

from adcircpy import AdcircMesh, AdcircRun
from adcircpy.forcing.waves.ww3 import WaveWatch3DataForcing
from adcircpy.forcing.winds.atmesh import AtmosphericMeshForcing
from adcircpy.server import SlurmConfig
from adcircpy.server.driver_file import DriverFile
from tests import download_mesh

MESH_URL = 'https://www.dropbox.com/s/1wk91r67cacf132/NetCDF_shinnecock_inlet.tar.bz2?dl=1'

DATA_DIRECTORY = Path(__file__).parent.absolute() / 'data'
REFERENCE_DIRECTORY = DATA_DIRECTORY / 'reference'
INPUT_DIRECTORY = DATA_DIRECTORY / 'NetCDF_Shinnecock_Inlet'
TEMPORARY_DIRECTORY = tempfile.TemporaryDirectory()
OUTPUT_DIRECTORY = Path(TEMPORARY_DIRECTORY.name)

INPUT_DIRECTORY.mkdir(exist_ok=True)

download_mesh(
    url=MESH_URL, directory=INPUT_DIRECTORY,
)


def test_slurm_driver():
    output_directory = OUTPUT_DIRECTORY / 'test_slurm_driver'
    reference_directory = REFERENCE_DIRECTORY / 'test_slurm_driver'
    output_directory.mkdir(parents=True, exist_ok=True)

    mesh = AdcircMesh.open(INPUT_DIRECTORY / 'fort.14', crs=4326)

    slurm = SlurmConfig(
        account='account',
        ntasks=1000,
        run_name='adcircpy/tests/test_configuration.py',
        partition='partition',
        walltime=timedelta(hours=8),
        mail_type='all',
        mail_user='example@email.gov',
        log_filename='test_configuration.log',
        modules=['intel/2020', 'impi/2020', 'netcdf/4.7.2-parallel'],
        path_prefix='$HOME/adcirc/build',
    )
    driver = AdcircRun(
        mesh=mesh,
        start_date=datetime.now(),
        end_date=timedelta(days=7),
        spinup_time=timedelta(days=5),
        server_config=slurm,
    )
    DriverFile(driver).write(output_directory / 'slurm.job', overwrite=True)

    with open(output_directory / 'slurm.job') as generated_file:
        with open(reference_directory / 'slurm.job') as reference_file:
            assert generated_file.read() == reference_file.read()


def test_configuration():
    output_directory = OUTPUT_DIRECTORY / 'test_configuration'
    reference_directory = REFERENCE_DIRECTORY / 'test_configuration'
    output_directory.mkdir(parents=True, exist_ok=True)

    mesh = AdcircMesh.open(INPUT_DIRECTORY / 'fort.14', crs=4326)

    spinup_time = timedelta(days=2)
    start_date = datetime(2015, 12, 14) + spinup_time
    end_date = start_date + timedelta(days=3)

    wind_forcing = AtmosphericMeshForcing(
        filename='Wind_HWRF_SANDY_Nov2018_ExtendedSmoothT.nc', nws=17, interval_seconds=3600,
    )
    wave_forcing = WaveWatch3DataForcing(
        filename='ww3.HWRF.NOV2018.2012_sxy.nc', nrs=5, interval_seconds=3600,
    )

    mesh.add_forcing(wind_forcing)
    mesh.add_forcing(wave_forcing)

    driver = AdcircRun(mesh, start_date, end_date, spinup_time,)

    driver.write(output_directory, overwrite=True)

    for reference_filename in reference_directory.iterdir():
        generated_filename = output_directory / reference_filename.name
        with open(generated_filename) as generated_file, open(
            reference_filename
        ) as reference_file:
            assert generated_file.readlines()[1:] == reference_file.readlines()[1:]
