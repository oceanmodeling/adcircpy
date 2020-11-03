from datetime import datetime, timedelta
import os
import pathlib
import tarfile
import unittest

import requests

from adcircpy import AdcircMesh, AdcircRun
from adcircpy.server import SlurmConfig
from adcircpy.server._driver_file import _DriverFile

DATA_DIRECTORY = pathlib.Path(__file__).parent.absolute() / 'data'
INPUT_DIRECTORY = DATA_DIRECTORY / 'input'
REFERENCE_DIRECTORY = DATA_DIRECTORY / 'reference'
OUTPUT_DIRECTORY = DATA_DIRECTORY / 'output'
FORT14_FILENAME = INPUT_DIRECTORY / 'fort.14'

os.makedirs(INPUT_DIRECTORY, exist_ok=True)
os.makedirs(OUTPUT_DIRECTORY, exist_ok=True)


class TestAdcircRun(unittest.TestCase):
    def setUp(self):
        # fetch Shinnecock Inlet test data
        if not FORT14_FILENAME.is_file():
            url = "https://www.dropbox.com/s/1wk91r67cacf132/" \
                  "NetCDF_shinnecock_inlet.tar.bz2?dl=1"
            remote_file = requests.get(url, stream=True)
            input_filename = DATA_DIRECTORY / 'NetCDF_shinnecock_inlet.tar.bz2'
            with open(input_filename, 'wb') as f:
                f.write(remote_file.raw.read())
            with tarfile.open(input_filename, 'r:bz2') as tar:
                tar.extractall(INPUT_DIRECTORY)
            os.remove(input_filename)

    def test_slurm_driver(self):
        output_directory = OUTPUT_DIRECTORY / 'slurm'
        reference_directory = REFERENCE_DIRECTORY / 'slurm'
        os.makedirs(output_directory, exist_ok=True)

        # open mesh file
        mesh = AdcircMesh.open(FORT14_FILENAME, crs=4326)

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
        _DriverFile(driver).write(output_directory / 'slurm.job',
                                  overwrite=True)

        with open(output_directory / 'slurm.job') as generated_file:
            with open(reference_directory / 'slurm.job') as reference_file:
                assert generated_file.read() == reference_file.read()


if __name__ == '__main__':
    unittest.main()
