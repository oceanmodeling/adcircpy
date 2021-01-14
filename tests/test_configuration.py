from datetime import datetime, timedelta
import os
import pathlib
import tarfile
import unittest

import numpy
import requests

from adcircpy import AdcircMesh, AdcircRun
from adcircpy.forcing.waves.ww3 import WaveWatch3DataForcing
from adcircpy.forcing.winds.atmesh import AtmosphericMeshForcing
from adcircpy.server import SlurmConfig
from adcircpy.server.driver_file import DriverFile

DATA_DIRECTORY = pathlib.Path(__file__).parent.absolute() / 'data'
REFERENCE_DIRECTORY = DATA_DIRECTORY / 'reference'
INPUT_DIRECTORY = DATA_DIRECTORY / 'input'
OUTPUT_DIRECTORY = DATA_DIRECTORY / 'output'
FORT14_FILENAME = INPUT_DIRECTORY / 'fort.14'

os.makedirs(INPUT_DIRECTORY, exist_ok=True)
os.makedirs(OUTPUT_DIRECTORY, exist_ok=True)


class TestAdcircRun(unittest.TestCase):
    def setUp(self):
        # fetch Shinnecock Inlet test data
        if not FORT14_FILENAME.is_file():
            url = 'https://www.dropbox.com/s/1wk91r67cacf132/NetCDF_shinnecock_inlet.tar.bz2?dl=1'
            remote_file = requests.get(url, stream=True)
            input_filename = DATA_DIRECTORY / 'NetCDF_shinnecock_inlet.tar.bz2'
            with open(input_filename, 'wb') as f:
                f.write(remote_file.raw.read())
            with tarfile.open(input_filename, 'r:bz2') as tar:
                tar.extractall(INPUT_DIRECTORY)
            os.remove(input_filename)

    def test_slurm_driver(self):
        output_directory = OUTPUT_DIRECTORY / 'test_slurm_driver'
        reference_directory = REFERENCE_DIRECTORY / 'test_slurm_driver'
        output_directory.mkdir(parents=True, exist_ok=True)

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
        DriverFile(driver).write(output_directory / 'slurm.job',
                                 overwrite=True)

        with open(output_directory / 'slurm.job') as generated_file:
            with open(reference_directory / 'slurm.job') as reference_file:
                assert generated_file.read() == reference_file.read()

    def test_configuration(self):
        output_directory = OUTPUT_DIRECTORY / 'test_configuration'
        reference_directory = REFERENCE_DIRECTORY / 'test_configuration'
        output_directory.mkdir(parents=True, exist_ok=True)

        # open mesh file
        mesh = AdcircMesh.open(FORT14_FILENAME, crs=4326)

        # let's generate the tau0 factor
        mesh.generate_tau0()

        # let's also add a mannings to the domain (constant for this example)
        mesh.mannings_n_at_sea_floor = numpy.full(mesh.values.shape, 0.025)

        # set simulation dates
        spinup_time = timedelta(days=2)
        start_date = datetime(2015, 12, 14) + spinup_time
        end_date = start_date + timedelta(days=3)

        wind_forcing = AtmosphericMeshForcing(17, 3600)
        wave_forcing = WaveWatch3DataForcing(5, 3600)

        mesh.add_forcing(wind_forcing)
        mesh.add_forcing(wave_forcing)

        # instantiate AdcircRun object.
        driver = AdcircRun(
                mesh,
                start_date,
                end_date,
                spinup_time,
        )

        driver.write(output_directory, overwrite=True)

        for reference_filename in reference_directory.iterdir():
            generated_filename = output_directory / reference_filename.name
            with open(generated_filename) as generated_file, \
                    open(reference_filename) as reference_file:
                assert generated_file.readlines()[1:] == \
                       reference_file.readlines()[1:]


if __name__ == '__main__':
    unittest.main()
