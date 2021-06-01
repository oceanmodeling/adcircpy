#!/usr/bin/env python
import pathlib
import tarfile
import tempfile
from unittest.mock import patch
import urllib.request

from adcircpy.cmd import tide_gen

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


def test_tide_gen():
    cmd = [
        'tide_gen',
        f'{FORT14.resolve()}',
        '2021-02-26T00:00:00',
        '15',
        '--mesh-crs=epsg:4326',
        '--output-file=/dev/null',
    ]

    with patch('sys.argv', cmd):
        tide_gen.main()
