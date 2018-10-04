#! /usr/bin/env python
import unittest
from AdcircPy.Datum import VDatum
from AdcircPyTests import AdcircPyEnvironment

class TestVDatum(AdcircPyEnvironment, unittest.TestCase):
  def setUp(self):
    super(TestVDatum, self).__init__()
    self.read_environment_variables()
  
  def test_VDatum_REST(self):
    pass


if __name__ == '__main__':
  unittest.main()

