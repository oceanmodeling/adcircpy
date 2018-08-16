from AdcircPy.core.Validation.COOPS import _HarmonicConstituents

products = ['air_gap',
             'air_pressure',
             'air_temperature',
             'conductivity',
             'currents',
             'currents_survey',
             'daily_mean',
             'datums',
             'high_low',
             'hourly_height',
             'humidity',
             'monthly_mean',
             'one_minute_water_level',
             'predictions',
             'salinity',
             'visibility',
             'water_level',
             'water_temperature',
             'wind']

class COOPS(dict):
  def __init__(self, **kwargs):
    self._products = products
    dict.__init__(self, **kwargs)



class HarmonicConstituents(COOPS):
  def __init__(self, **kwargs):
    COOPS.__init__(self, **kwargs)
  
  @staticmethod
  def from_stations_list(stations):
    return _HarmonicConstituents._from_stations_list(stations)
