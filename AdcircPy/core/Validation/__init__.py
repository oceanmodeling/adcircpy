
class USGS(dict):
  def __init__(self, **kwargs):
    self.epsg   = kwargs.pop("epsg", None)
    dict.__init__(self, **kwargs)

class COOPS(dict):
  def __init__(self, **kwargs):
    self._products = ['air_gap',
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
    dict.__init__(self, **kwargs)