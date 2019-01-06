#! /usr/bin/env python
from AdcircPyTests import AdcircPyEnvironment
import unittest
from datetime import datetime, timedelta
from AdcircPy.Tides import TidalForcing

class TestTidalForcing(unittest.TestCase):
  """
  This class initializes tidal forcings information required to
  run some tidal models.
  Note:
    TPXO initializations is not part of TidalForcing because
    it depends on the mesh boundaries, therefore TPXO initialization
    is called on the _AdcircRun class.
  """

  def setUp(self):
    self.now = datetime.now()
    self.later = self.now + timedelta(days=7) 

  def test_TidalForcing(self):
    TidalForcing(self.now, self.later)

  def test_wrong_spinup_datetime_raises(self):
    wrong_spinup_date = self.now + timedelta(days=5)
    with self.assertRaises(Exception) as context:
      TidalForcing(self.now, self.later, spinup_date=wrong_spinup_date)
    self.assertTrue('spinup_date must be smaller than start_date.' in str(context.exception))

  def test_wrong_spinup_date_datatype_raise(self):
    wrong_spinup_date_datatype = ''
    with self.assertRaises(Exception) as context:
      TidalForcing(self.now, self.later, spinup_date=wrong_spinup_date_datatype)
    self.assertTrue('spinup_date must be a datetime.datetime instance.' in str(context.exception))

  def test_wrong_end_date_datatype_raise(self):
    wrong_end_date_datatype = ''
    with self.assertRaises(Exception) as context:
      TidalForcing(self.now, wrong_end_date_datatype)
    self.assertTrue('end_date must be a datetime.datetime instance.' in str(context.exception))


  def test_wrong_start_date_datatype_raise(self):
    wrong_start_date_datatype = ''
    with self.assertRaises(Exception) as context:
      TidalForcing(wrong_start_date_datatype, self.later)
    self.assertTrue('start_date must be a datetime.datetime instance.' in str(context.exception))

if __name__ == '__main__':
  unittest.main()