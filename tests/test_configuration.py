# ! /usr/bin/env python
from copy import copy
from datetime import datetime, timedelta

from adcircpy import AdcircMesh, AdcircRun
from adcircpy.forcing.waves.ww3 import WaveWatch3DataForcing
from adcircpy.forcing.winds.atmesh import AtmosphericMeshForcing
from adcircpy.fort15 import StationType
from adcircpy.server import SlurmConfig
from adcircpy.server.driver_file import DriverFile

# noinspection PyUnresolvedReferences
from tests import (
    check_reference_directory,
    INPUT_DIRECTORY,
    OUTPUT_DIRECTORY,
    REFERENCE_DIRECTORY,
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


def test_import_stations(shinnecock_mesh_directory):
    input_directory = INPUT_DIRECTORY / 'test_import_stations'

    mesh = AdcircMesh.open(shinnecock_mesh_directory / 'fort.14', crs=4326)

    spinup_time = timedelta(days=2)
    start_date = datetime(2015, 12, 14) + spinup_time
    end_date = start_date + timedelta(days=3)

    driver_1 = AdcircRun(copy(mesh), start_date, end_date, spinup_time)
    driver_1.import_stations(input_directory / 'stations_1.txt')

    driver_2 = AdcircRun(copy(mesh), start_date, end_date, spinup_time)
    driver_2.import_stations(input_directory / 'stations_2.txt')

    driver_3 = AdcircRun(copy(mesh), start_date, end_date, spinup_time)
    driver_3.import_stations(input_directory / 'stations_3.txt')

    driver_4 = AdcircRun(copy(mesh), start_date, end_date, spinup_time)
    driver_4.import_stations(
        input_directory / 'stations_3.txt',
        station_types=['elevation', 'NSTAC', StationType.METEOROLOGICAL],
    )

    assert driver_1.elevation_stations == {'8512769': (-72.5772, 40.823)}
    assert driver_1.velocity_stations == {}
    assert driver_1.concentration_stations == {}
    assert driver_1.meteorological_stations == {}

    assert driver_2.elevation_stations == {'8512769': (-72.5772, 40.823)}
    assert driver_2.velocity_stations == {'8512769': (-72.5772, 40.823)}
    assert driver_2.concentration_stations == {}
    assert driver_2.meteorological_stations == {}

    assert driver_3.elevation_stations == {'8512769': (-72.5772, 40.823)}
    assert driver_3.velocity_stations == {'8512769': (-72.5772, 40.823)}
    assert driver_3.concentration_stations == {}
    assert driver_3.meteorological_stations == {}

    assert driver_4.elevation_stations == {'8512769': (-72.5772, 40.823)}
    assert driver_4.velocity_stations == {}
    assert driver_4.concentration_stations == {'8512769': (-72.5772, 40.823)}
    assert driver_4.meteorological_stations == {}
