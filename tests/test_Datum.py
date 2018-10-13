#! /usr/bin/env python
import unittest
from AdcircPy.Datum import DatumMesh
from AdcircPyTests import AdcircPyEnvironment

class TestDatum(AdcircPyEnvironment, unittest.TestCase):
  def setUp(self):
    super(TestDatum, self).__init__()
    self.read_environment_variables()
  
  def test_generate_DatumMesh(self):
    grid = DatumGrid.build_datum_mesh(self.os.getenv('FORT14_PATH'),
                                      source_hdatum='NAD83', target_hdatum='NAD83',
                                      source_vdatum="LMSL", target_vdatum="MHHW",
                                      vdatum_jar_path=self._os.getenv('VDATUM_JAR_PATH'))
    grid.dump('./fort14_lmsl2mhhw.grd')
    self.os.remove('./fort14_lmsl2mhhw.grd')

if __name__ == '__main__':
  unittest.main()
