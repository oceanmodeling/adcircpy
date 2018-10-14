#! /usr/bin/env python
import unittest
from datetime import datetime, timedelta
from AdcircPy.Winds import BestTrack

class TestBestTrack(unittest.TestCase):
  """
  unittest class for BestTrack wind generator.
  """
  def setUp(self):
    self.Sandy_id = 'AL182012'

  def test_BestTrack_Sandy_full(self):
    bt = BestTrack(self.Sandy_id)
    bt.printf()

  def test_BestTrack_Sandy_official(self):
    start_date = datetime(2012, 10, 26)
    end_date = datetime(2012, 10, 31, 12)
    bt = BestTrack(self.Sandy_id, start_date, end_date)
    bt.printf()

if __name__ == '__main__':
  unittest.main()