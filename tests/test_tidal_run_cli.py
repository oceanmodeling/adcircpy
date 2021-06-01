#!/usr/bin/env python
from datetime import datetime
import pathlib
import shutil
import tarfile
import tempfile
from unittest.mock import patch
import urllib.request

from adcircpy.cmd import tidal_run

DATA_DIRECTORY = pathlib.Path(__file__).parent.absolute() / 'data'
FORT14 = DATA_DIRECTORY / 'NetCDF_Shinnecock_Inlet/fort.14'

if not FORT14.is_file():
    url = 'https://www.dropbox.com/s/1wk91r67cacf132/'
    url += 'NetCDF_shinnecock_inlet.tar.bz2?dl=1'
    g = urllib.request.urlopen(url)
    tmpfile = tempfile.NamedTemporaryFile()
    with open(tmpfile.name, 'b+w') as f:
        f.write(g.read())
    with tarfile.open(tmpfile.name, 'r:bz2') as tar:
        tar.extractall(DATA_DIRECTORY / 'NetCDF_Shinnecock_Inlet')


def test_tidal_run():
    cmd = [
        'tidal_run',
        f'{FORT14.resolve()}',
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
        cmd.append('--skip-run')

    with patch('sys.argv', cmd):
        tidal_run.main()
