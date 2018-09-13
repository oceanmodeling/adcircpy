
class _OutputStations(dict):
  def __init__(self, time, **station_data):
    self.time = time
    dict.__init__(self, **station_data)

class ElevationStations(_OutputStations):
  def __init__(self, time, **station_data):
    _OutputStations.__init__(self, time, **station_data)
  
  @staticmethod
  def from_netcdf(path):
    return _ElevationStations._from_netcdf(_OutputStations)

class HarmonicConstituentsStations(dict):
  def __init__(self, **kwargs):
    dict.__init__(self, **kwargs)

  @staticmethod
  def from_fort51(path, fort14, fort15):
    return _HarmonicConstituentsStations._from_fort51(path, fort14, fort15)


class VelocityStations(_OutputStations):
  def __init__(self, **kwargs):
    _OutputStations.__init__(self, **kwargs)

  def make_plot(self, station, **kwargs):
    return _VelocityStations._make_plot(self, station, **kwargs)

  @staticmethod
  def from_netcdf(path):
    return _VelocityStations._from_netcdf(path)
