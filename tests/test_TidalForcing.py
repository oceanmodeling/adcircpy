#!/usr/bin/env python
import unittest
from datetime import datetime, timedelta
from adcircpy.model import TidalForcing


class TidalForcingTestCase(unittest.TestCase):

    def test_control_Sandy(self):
        """
        K1    0.94843     298.0830551604721
        O1    0.91613     167.70941789973676
        P1    1.0         70.01137799920002
        Q1    0.91613     266.0233347217413
        N2    1.0203      208.1293472169341
        M2    1.0203      109.81543039492954
        S2    1.0         0.0
        K2    0.86465     55.559822189650276
        M4    1.041       219.63086078985907
        MS4   1.0203      109.81543039492954
        """
        t = TidalForcing()
        t.start_date = datetime(2012, 10, 26, 00)
        t.end_date = t.start_date + timedelta(days=4.25)
        t.spinup_time = timedelta(days=15)
        ordered = ["K1", "O1", "P1", "Q1", "N2", "M2", "S2", "K2",
                   "Mf", "Mm", "M4", "MS4"]
        for constituent in ordered:
            t.use_constituent(constituent)
        for constituent in ordered:
            s = '{:<6}'.format(constituent)
            s += '{:<7.5}'.format(t.get_nodal_factor(constituent))
            s += 5*' '
            s += '{} '.format(t.get_greenwich_term(constituent))
            print(s)

    def _test_Isabel2003(self):
        tidal_forcing = TidalForcing()
        tidal_forcing.start_date = datetime(2003, 9, 7, 12)
        tidal_forcing.end_date = datetime(2003, 9, 19, 0)
        tidal_forcing.spinup_time = timedelta(days=15)
        for forcing in tidal_forcing:
            print(forcing)


if __name__ == '__main__':
    unittest.main()
