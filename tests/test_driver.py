from datetime import datetime, timedelta
import os
import pathlib
import tarfile
import unittest

import requests

from adcircpy import AdcircMesh, AdcircRun, Tides
from adcircpy.server import SlurmConfig

DATA_DIRECTORY = pathlib.Path(__file__).parent.absolute() / 'data'
INPUT_DIRECTORY = DATA_DIRECTORY / 'input'
OUTPUT_DIRECTORY = DATA_DIRECTORY / 'output'
REFERENCE_DIRECTORY = DATA_DIRECTORY / 'reference'
FORT14_FILENAME = INPUT_DIRECTORY / 'fort.14'

os.makedirs(DATA_DIRECTORY, exist_ok=True)
os.makedirs(INPUT_DIRECTORY, exist_ok=True)
os.makedirs(OUTPUT_DIRECTORY, exist_ok=True)


class TestAdcircRun(unittest.TestCase):
    def test_slurm_driver(self):
        # fetch Shinnecock Inlet test data
        if not FORT14_FILENAME.is_file():
            url = "https://www.dropbox.com/s/1wk91r67cacf132/" \
                  "NetCDF_shinnecock_inlet.tar.bz2?dl=1"
            remote_file = requests.get(url, stream=True)
            temp_filename = DATA_DIRECTORY / 'temp.tar.gz'
            with open(temp_filename, 'wb') as f:
                f.write(remote_file.raw.read())
            with tarfile.open(temp_filename, 'r:bz2') as tar:
                tar.extractall(INPUT_DIRECTORY)
            os.remove(temp_filename)

        # open mesh file
        mesh = AdcircMesh.open(FORT14_FILENAME, crs=4326)

        # init tidal forcing and setup requests
        tidal_forcing = Tides()
        tidal_forcing.use_all()

        mesh.add_forcing(tidal_forcing)

        # instantiate AdcircRun object.
        slurm = SlurmConfig(
            account='account',
            ntasks=1000,
            run_name='AdcircPy/examples/example_3.py',
            partition='partition',
            walltime=timedelta(hours=8),
            mail_type='all',
            mail_user='example@email.gov',
            log_filename='example_3.log',
            modules=['intel/2020', 'impi/2020', 'netcdf/4.7.2-parallel'],
            path_prefix='$HOME/adcirc/build'
        )
        driver = AdcircRun(
            mesh=mesh,
            start_date=datetime.now(),
            end_date=timedelta(days=7),
            spinup_time=timedelta(days=5),
            server_config=slurm
        )
        driver.write(OUTPUT_DIRECTORY / 'slurm', overwrite=True)

        with open(OUTPUT_DIRECTORY / 'slurm/slurm.job') as generated:
            with open(REFERENCE_DIRECTORY / 'slurm/slurm.job') as reference:
                assert generated.read() == reference.read()


if __name__ == '__main__':
    unittest.main()
