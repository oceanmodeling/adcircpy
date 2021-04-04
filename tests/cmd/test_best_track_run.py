#!/usr/bin/env python
import pathlib
import shutil
import tarfile
import tempfile
import unittest
import urllib.request

from adcircpy.cmd.argument_parser import get_parser
from adcircpy.cmd.best_track_run import BestTrackRunCommand


DATA_DIRECTORY = pathlib.Path(__file__).parent.parent.absolute() / 'data'
FORT14 = DATA_DIRECTORY / "NetCDF_Shinnecock_Inlet/fort.14"


class BestTrackRunCliTestCase(unittest.TestCase):

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

    def test_best_track_run(self):
        cmd = [
                f'{FORT14.resolve()}',
                'Sandy2012',
                '--spinup-days=0.5',
                '--crs=EPSG:4326',
                '--output-directory=/tmp/test',
                '--constituents=all',
                '--overwrite',
                '--timestep=10.',
                '--tau0-gen',
                f'--stations-file={pathlib.Path(__file__).parent}/stations.txt',
                '--elev-stat=6.',
                '--run-days=0.5',
                # '--nproc=1'
            ]
        if shutil.which('padcirc') is None:
            cmd.append('--skip-run')

        BestTrackRunCommand(
            get_parser('best_track').parse_args(cmd)).run()


if __name__ == '__main__':
    unittest.main()
