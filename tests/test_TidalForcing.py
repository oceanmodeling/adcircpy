#!/usr/bin/env python
import unittest
from datetime import datetime, timedelta
from adcircpy.model import TidalForcing


class TidalForcingTestCase(unittest.TestCase):

    def test_Isabel2003(self):
        tidal_forcing = TidalForcing()
        tidal_forcing.start_date = datetime(2003, 9, 7, 12)
        tidal_forcing.end_date = datetime(2003, 9, 19, 0)
        tidal_forcing.spinup_time = timedelta(days=15)
        for forcing in tidal_forcing:
            print(forcing)


if __name__ == '__main__':
    unittest.main()
