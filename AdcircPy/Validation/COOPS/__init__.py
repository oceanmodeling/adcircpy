from AdcircPy.Validation.COOPS import _TidalStations
from AdcircPy.Validation.COOPS import _HarmonicConstituents

class _REST(dict):
  
  _params={"format"       : "json",
           "units"         : "metric", 
           "time_zone"    : "gmt",
           "application" : "AdcircPy",
           "datum"       : "msl"}
  
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
    self._format = 'json'
    self._unit = 'metric'
    self._time_zone = 'gmt'
    self._application = 'AdcircPy'
    self._url = "https://tidesandcurrents.noaa.gov/api/datagetter?"
    
class TidalStations(_REST):
  def __init__(self, **kwargs):
    dict.__init__(self, **kwargs)

  @staticmethod
  def from_station_list(stations, start_date, end_date):
    return _TidalStations._from_station_list(stations, start_date, end_date)

class HarmonicConstituents(dict):
  def __init__(self, **kwargs):
    dict.__init__(self, **kwargs)
  
  @staticmethod
  def from_station_list(stations):
    return _HarmonicConstituents._from_station_list(stations)

