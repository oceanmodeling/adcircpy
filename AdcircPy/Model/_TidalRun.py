# local imports
from AdcircPy.Model._AdcircRun import _AdcircRun
from AdcircPy.Model._TidalForcing import _TidalForcing


class _TidalRun(_AdcircRun):
    def __init__(self, AdcircMesh, start_date, end_date, spinup_days=0.,
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
