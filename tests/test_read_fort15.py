#! /usr/bin/env python
import unittest
from AdcircPy import AdcircPy
from AdcircPyTests import AdcircPyEnvironment

class testFort15(AdcircPyEnvironment, unittest.TestCase):
  def setUp(self):
    super(testFort15, self).__init__()
    self.read_environment_variables()
  
  def test_read_coldstart(self):
    AdcircPy.read_mesh(fort14=self._os.getenv("FORT14_PATH"),
                       fort15=self._os.getenv("FORT15_COLDSTART_PATH"))

  def test_read_hotstart(self):
    AdcircPy.read_mesh(fort14=self._os.getenv("FORT14_PATH"),
                       fort15=self._os.getenv("FORT15_HOTSTART_PATH"))

if __name__ == '__main__':
  unittest.main()

