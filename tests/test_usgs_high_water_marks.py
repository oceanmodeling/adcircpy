#! /usr/bin/env python
import unittest
from AdcircPyTests import AdcircPyEnvironment
from AdcircPy.Validation import HighWaterMarks

class testUSGSHighWaterMarks(AdcircPyEnvironment, unittest.TestCase):
  def setUp(self):
    super(testUSGSHighWaterMarks, self).__init__()
    self.read_environment_variables()

  def test_usgs_high_water_mark_import(self):
    filter_dict={"poor":True, "fair":True, "riverine":True, "non_still_water":True}
    target_datum='LMSL'
    HWM = HighWaterMarks.from_event_name(self._os.getenv("STORM_EVENT_NAME"), target_datum, filter_dict=filter_dict)

if __name__ == '__main__':
  unittest.main()
