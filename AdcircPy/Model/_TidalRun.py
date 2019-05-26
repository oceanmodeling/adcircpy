# local imports
from AdcircPy.Model._AdcircRun import _AdcircRun
from AdcircPy.Model._TidalForcing import _TidalForcing

# unittest imports
import unittest
import os
from datetime import datetime, timedelta
from AdcircPy import Model


class _TidalRun(_AdcircRun):
    def __init__(self, AdcircMesh, start_date, end_date, spinup_days=7.,
                 constituents='all', netcdf=True, ElevationGlobalOutput=None,
                 VelocityGlobalOutput=None, MeteorologicalGlobalOutput=None,
                 ElevationStationsOutput=None, VelocityStationsOutput=None,
                 MeteorologicalStationsOutput=None, **fort15):
        TidalForcing = _TidalForcing(start_date, end_date, spinup_days)
        super(_TidalRun, self).__init__(AdcircMesh, TidalForcing, None, netcdf,
                                        ElevationGlobalOutput,
                                        VelocityGlobalOutput,
                                        MeteorologicalGlobalOutput,
                                        ElevationStationsOutput,
                                        VelocityStationsOutput,
                                        MeteorologicalStationsOutput,
                                        **fort15)


class TidalRunTestCase(unittest.TestCase):
    def setUp(self):
        self.AdcircMesh = Model.AdcircMesh(os.getenv('ShinnecockInlet'), 4326,
                                           'LMSL')
        self.ElevationStationsOutput = Model.ElevationStationsOutput \
            .from_fort15(os.getenv('FORT15_HOTSTART_PATH'),
                         sampling_frequency=timedelta(minutes=6),
                         harmonic_analysis=True)
        self.ElevationGlobalOutput = Model.ElevationGlobalOutput(
                                    sampling_frequency=timedelta(minutes=15),
                                    harmonic_analysis=True)

    def test_TidalRun_90dayTidal(self):
        # Datetime is the same as the HSOFS offical 90 day tidal run.
        spinup_date = datetime(2013, 8, 1, 0, 0, 0)
        start_time = spinup_date + timedelta(days=30)
        end_time = spinup_date + timedelta(days=90)
        TidalRun = self.AdcircMesh.TidalRun(
                        start_time, end_time, spinup_days=7.,
                        ElevationStationsOutput=self.ElevationStationsOutput,
                        ElevationGlobalOutput=self.ElevationGlobalOutput)
        TidalRun.dump()


if __name__ == '__main__':
    unittest.main()
