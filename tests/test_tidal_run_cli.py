#! /usr/bin/env python

from datetime import datetime
from pathlib import Path
import shutil
from unittest.mock import patch

from adcircpy.cmd import tidal_run
from tests import download_mesh

DATA_DIRECTORY = Path(__file__).parent.absolute().resolve() / 'data'
INPUT_DIRECTORY = DATA_DIRECTORY / 'input' / 'NetCDF_Shinnecock_Inlet'

download_mesh(
    url='https://www.dropbox.com/s/1wk91r67cacf132/NetCDF_shinnecock_inlet.tar.bz2?dl=1',
    directory=INPUT_DIRECTORY,
)


def test_tidal_run():
    cmd = [
        'tidal_run',
        f'{INPUT_DIRECTORY / "fort.14"}',
        f"{datetime.strftime(datetime.utcnow(), '%Y-%m-%dT%H:%M:%S')}",
        '0.5',
        '--spinup-days=0.5',
        '--crs=EPSG:4326',
        '--output-directory=/tmp/test',
        '--constituents=all',
        '--overwrite',
        '--timestep=10.',
    ]
    if shutil.which('padcirc') is None:
        if shutil.which('adcirc') is None:
            cmd.append('--nproc=1')
        else:
            cmd.append('--skip-run')

    with patch('sys.argv', cmd):
        tidal_run.main()
