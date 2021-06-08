#! /usr/bin/env python

from adcircpy.cmd import tide_gen

# noinspection PyUnresolvedReferences
from tests import (
    OUTPUT_DIRECTORY,
    REFERENCE_DIRECTORY,
    check_reference_directory,
    shinnecock_mesh_directory,
)


def test_tide_gen(shinnecock_mesh_directory, mocker):
    output_directory = OUTPUT_DIRECTORY / 'test_tide_gen'
    reference_directory = REFERENCE_DIRECTORY / 'test_tide_gen'

    if not output_directory.exists():
        output_directory.mkdir(parents=True, exist_ok=True)

    cmd = [
        'tide_gen',
        f'{shinnecock_mesh_directory / "fort.14"}',
        '2021-02-26T00:00:00',
        '15',
        '--mesh-crs=epsg:4326',
        f'--output-file={output_directory / "fort.15"}',
    ]
    mocker.patch('sys.argv', cmd)

    tide_gen.main()

    check_reference_directory(
        output_directory, reference_directory, skip_lines={'fort.15': [0, -1]}
    )
