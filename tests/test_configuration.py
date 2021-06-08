# ! /usr/bin/env python

from datetime import datetime, timedelta

from adcircpy import AdcircMesh, AdcircRun
from adcircpy.forcing.waves.ww3 import WaveWatch3DataForcing
from adcircpy.forcing.winds.atmesh import AtmosphericMeshForcing
from adcircpy.server import SlurmConfig
from adcircpy.server.driver_file import DriverFile

# noinspection PyUnresolvedReferences
from tests import (
    OUTPUT_DIRECTORY,
    REFERENCE_DIRECTORY,
    check_reference_directory,
    shinnecock_mesh_directory,
)


def test_slurm_driver(shinnecock_mesh_directory):
    output_directory = OUTPUT_DIRECTORY / 'test_slurm_driver'
    reference_directory = REFERENCE_DIRECTORY / 'test_slurm_driver'

    if not output_directory.exists():
        output_directory.mkdir(parents=True, exist_ok=True)

    mesh = AdcircMesh.open(shinnecock_mesh_directory / 'fort.14', crs=4326)

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

    check_reference_directory(output_directory, reference_directory)


def test_configuration(shinnecock_mesh_directory):
    output_directory = OUTPUT_DIRECTORY / 'test_configuration'
    reference_directory = REFERENCE_DIRECTORY / 'test_configuration'

    if not output_directory.exists():
        output_directory.mkdir(parents=True, exist_ok=True)

    mesh = AdcircMesh.open(shinnecock_mesh_directory / 'fort.14', crs=4326)

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

    driver.write(output_directory, overwrite=True, nproc=2)

    check_reference_directory(
        output_directory, reference_directory, skip_lines={'fort.15': [0]},
    )
