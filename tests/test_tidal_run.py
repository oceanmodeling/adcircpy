#! /usr/bin/env python

from datetime import datetime, timedelta
import shutil

from adcircpy.cmd import tidal_run
from adcircpy.driver import AdcircRun
from adcircpy.forcing import Tides
from adcircpy.mesh import AdcircMesh

# noinspection PyUnresolvedReferences
from tests import (
    OUTPUT_DIRECTORY,
    REFERENCE_DIRECTORY,
    check_reference_directory,
    shinnecock_mesh_directory,
)


def test_tidal_run(shinnecock_mesh_directory):
    output_directory = OUTPUT_DIRECTORY / 'test_tidal_run'
    reference_directory = REFERENCE_DIRECTORY / 'test_tidal_run'

    if not output_directory.exists():
        output_directory.mkdir(parents=True, exist_ok=True)

    mesh = AdcircMesh.open(shinnecock_mesh_directory / 'fort.14', crs=4326)

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

    if shutil.which('padcirc') is not None:
        driver.run(output_directory, nproc=2, overwrite=True)
    else:
        driver.write(output_directory, nproc=2, overwrite=True)

    check_reference_directory(
        output_directory,
        reference_directory,
        skip_lines={'fort.15': [0, *range(32, 47, 2), *range(49, 64, 2), -2]},
    )


def test_tidal_run_cli(shinnecock_mesh_directory, mocker):
    output_directory = OUTPUT_DIRECTORY / 'test_tidal_run_cli'
    reference_directory = REFERENCE_DIRECTORY / 'test_tidal_run_cli'

    if not output_directory.exists():
        output_directory.mkdir(parents=True, exist_ok=True)

    cmd = [
        'tidal_run',
        f'{shinnecock_mesh_directory / "fort.14"}',
        f"{datetime.strftime(datetime.utcnow(), '%Y-%m-%dT%H:%M:%S')}",
        '0.5',
        '--spinup-days=0.5',
        '--crs=EPSG:4326',
        f'--output-directory={output_directory}',
        '--constituents=all',
        '--overwrite',
        '--timestep=10.',
        '--nproc=2',
    ]
    if shutil.which('padcirc') is None:
        if shutil.which('adcirc') is not None:
            cmd.append('--nproc=1')
        else:
            cmd.append('--skip-run')
    mocker.patch('sys.argv', cmd)

    tidal_run.main()

    check_reference_directory(
        output_directory,
        reference_directory,
        skip_lines={'fort.15': [0, *range(32, 47, 2), *range(49, 64, 2), -2]},
    )
