#!/usr/bin/env python
from datetime import datetime, timedelta
import pathlib
import shutil
import tarfile
import tempfile
import unittest
import urllib.request

from adcircpy.mesh import AdcircMesh
from adcircpy.driver import AdcircRun
from adcircpy.forcing import Tides


DATA_DIRECTORY = pathlib.Path(__file__).parent.absolute() / 'data'
FORT14 = DATA_DIRECTORY / "NetCDF_Shinnecock_Inlet/fort.14"


class TidalRunTestCase(unittest.TestCase):

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
        # open mesh file
        mesh = AdcircMesh.open(FORT14, crs=4326)

        # init tidal forcing and setup requests
        tidal_forcing = Tides()
        tidal_forcing.use_all()

        mesh.add_forcing(tidal_forcing)
        now = datetime.utcnow()
        driver = AdcircRun(
            mesh,
            start_date=now,
            end_date=now+timedelta(days=0.5),
            spinup_time=timedelta(days=0.5),
        )
        driver.timestep = 10.
        if shutil.which('padcirc') is None:
            tmpdir = tempfile.TemporaryDirectory()
            driver.write(tmpdir.name)
        else:
            driver.run()


if __name__ == '__main__':
    unittest.main()
