#!/usr/bin/env python
from datetime import datetime
import pathlib
import sys
import shutil
import tarfile
import tempfile
import unittest
from unittest.mock import patch
import urllib.request

from adcircpy.cmd import tidal_run


DATA_DIRECTORY = pathlib.Path(__file__).parent.absolute() / 'data'
FORT14 = DATA_DIRECTORY / "NetCDF_Shinnecock_Inlet/fort.14"


class TidalRunCliTestCase(unittest.TestCase):

    def setUp(self):
        if not FORT14.is_file():
            url = "https://www.dropbox.com/s/1wk91r67cacf132/"
            url += "NetCDF_shinnecock_inlet.tar.bz2?dl=1"
            g = urllib.request.urlopen(url)
            tmpfile = tempfile.NamedTemporaryFile()
            with open(tmpfile.name, 'b+w') as f:
                f.write(g.read())
            with tarfile.open(tmpfile.name, "r:bz2") as tar:
                tar.extractall(DATA_DIRECTORY / "NetCDF_Shinnecock_Inlet")

    def test_tidal_run(self):
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
        with patch.object(sys, 'argv', cmd):
            tidal_run.main()


if __name__ == '__main__':
    unittest.main()
