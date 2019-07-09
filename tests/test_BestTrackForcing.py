#!/usr/bin/env python
import unittest
from datetime import datetime
from adcircpy.model import BestTrackForcing


class BestTrackForcingTestCase(unittest.TestCase):

    def test_Isabel2003(self):
        #  ----------- add wind forcing
        WindForcing = BestTrackForcing()
        WindForcing.storm_id = 'AL132003'
        WindForcing.remove_TS()
        WindForcing.remove_EX()
        print(WindForcing.fort22)

    def _test_Sandy2012(self):
        #  ----------- add wind forcing
        WindForcing = BestTrackForcing()
        WindForcing.storm_id = 'AL182012'
        WindForcing.start_date = datetime(2012, 10, 27, 12)
        WindForcing.end_date = datetime(2012, 10, 29, 12)
        print(WindForcing.fort22)
        # WindForcing.remove_TS()
        # WindForcing.remove_EX()


if __name__ == '__main__':
    unittest.main()
