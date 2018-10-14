#! /usr/bin/env python
import unittest
from datetime import datetime, timedelta
from AdcircPy.Winds import BestTrack

class TestBestTrack(unittest.TestCase):
  """
  unittest class for BestTrack wind generator.
  """

  def test_BestTrack_Sandy(self):
    hurdat2_id = 'AL182012'
    bt = BestTrack(hurdat2_id)
    bt.printf()


if __name__ == '__main__':
  unittest.main()