#! /usr/bin/env python
import unittest
from AdcircPyTests import AdcircPyEnvironment
from datetime import datetime, timedelta
from AdcircPy import AdcircPy
from AdcircPy import ElevationStationsOutput as ESO
from AdcircPy import ElevationGlobalOutput as EGO

class GenerateSandyHindcast(AdcircPyEnvironment, unittest.TestCase):
  def setUp(self):
    super(GenerateSandyHindcast, self).__init__()
    self.read_environment_variables()
    self.AdcircMesh = AdcircPy.read_mesh(fort14=self.os.getenv('FORT14_PATH'), fort13=self.os.getenv('FORT13_PATH'))
    self.ElevationStationsOutput = ESO.from_fort15(self.os.getenv('FORT15_HOTSTART_PATH'))
    self.ElevationGlobalOutput = EGO()
    self.hurdat2_id = 'AL182012'

  def test_generate_90DayTidal(self):
    # Datetime is the same as the HSOFS offical Sandy hindcast run.
    HindcastRun = self.AdcircMesh.HindcastRun(self.hurdat2_id,
                                          ElevationStationsOutput=self.ElevationStationsOutput,
                                          ElevationGlobalOutput=self.ElevationGlobalOutput)
    HindcastRun.dump("./", printf=True)
    self.os.remove('./fort.15.coldstart')
    self.os.remove('./fort.15.hotstart')

if __name__ == '__main__':
    unittest.main()