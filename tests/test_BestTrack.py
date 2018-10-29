#! /usr/bin/env python
import unittest
from datetime import datetime, timedelta
from AdcircPy.Winds import BestTrack

class TestBestTrack(unittest.TestCase):
  """
  unittest class for BestTrack wind generator.
  """
  def test_BestTrack_Sandy_full(self):
    bt = BestTrack('AL182012')
    bt.printf()
    bt.plot_track()
    bt.plt.show()

  def test_BestTrack_Sandy_official(self):
    start_date = datetime(2012, 10, 26, 0)
    end_date = datetime(2012, 10, 31, 6)
    bt = BestTrack('AL182012', start_date, end_date)
    bt.printf()
    bt.plot_track()
    bt.plt.show()

  def test_BestTrack_Charley_full(self):
    bt = BestTrack('AL032004')
    bt.printf()
    bt.plot_track()
    bt.plt.show()



if __name__ == '__main__':
  unittest.main()