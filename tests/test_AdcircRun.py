import unittest
import os
import shutil
from datetime import datetime, timedelta
from AdcircPy import AdcircPy
from AdcircPy import TidalForcing
from AdcircPy import AdcircRun
from AdcircPy import ElevationStationsOutput
from AdcircPyTests import AdcircPyEnvironment


class GenerateAdcircRunTests(AdcircPyEnvironment, unittest.TestCase):
  def setUp(self):
    super(GenerateAdcircRunTests, self).__init__()
    self.read_environment_variables()
    self.AdcircMesh = AdcircPy.read_mesh(fort14=self._os.getenv('FORT14_PATH'),
                                         # fort13=self._os.getenv('FORT13_PATH')
                                         )
  
  def test_generate_90DayTidal(self):
    # Datetime is the same as the HSOFS offical 90 day tidal run.
    spinup_date  = datetime(2013, 8, 1, 0, 0, 0)
    start_time   = spinup_date + timedelta(days=30)
    end_time     = spinup_date + timedelta(days=90)
    # elevationStationOutputs = ElevationStationsOutput.from_fort15(self._os.getenv('FORT15_HOTSTART_PATH'))
    tidalForcing = TidalForcing(start_time, end_time,
                                constituents=None,
                                spinup_date=spinup_date)
    modelRun     = AdcircRun(self.AdcircMesh, Tides=tidalForcing)
    modelRun.dump("./")
    with open('./fort.15.coldstart', 'r') as fort15:
        for line in fort15.read().splitlines():
            print(line)

    with open('./fort.15.hotstart', 'r') as fort15:
        for line in fort15.read().splitlines():
            print(line)
    os.remove('./fort.15.coldstart')
    os.remove('./fort.15.hotstart')

if __name__ == '__main__':
    unittest.main()