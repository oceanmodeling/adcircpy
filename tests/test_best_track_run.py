#! /usr/bin/env python

import shutil

from adcircpy.cmd import best_track_run

# noinspection PyUnresolvedReferences
from tests import (
    check_reference_directory,
    INPUT_DIRECTORY,
    OUTPUT_DIRECTORY,
    REFERENCE_DIRECTORY,
    shinnecock_mesh_directory,
)


def test_best_track_run(shinnecock_mesh_directory, mocker):
    input_directory = INPUT_DIRECTORY / 'test_best_track_run'
    output_directory = OUTPUT_DIRECTORY / 'test_best_track_run'
    reference_directory = REFERENCE_DIRECTORY / 'test_best_track_run'

    if not output_directory.exists():
        output_directory.mkdir(parents=True, exist_ok=True)

    cmd = [
        'best_track_run',
        f'{shinnecock_mesh_directory / "fort.14"}',
        'Sandy2012',
        '--spinup-days=0.5',
        '--crs=EPSG:4326',
        f'--output-directory={str(output_directory)}',
        '--constituents=all',
        '--overwrite',
        '--timestep=10.',
        '--tau0-gen',
        f'--stations-file={input_directory / "stations.txt"}',
        '--elev-stat=6.',
        '--run-days=0.5',
        '--nproc=2',
    ]
    if shutil.which('padcirc') is None:
        if shutil.which('adcirc') is not None:
            cmd.append('--nproc=1')
        else:
            cmd.append('--skip-run')
    mocker.patch('sys.argv', cmd)

    best_track_run.main()

    check_reference_directory(
        output_directory, reference_directory, skip_lines={'fort.15': [0, -1]}
    )
