class VelocityStations(_OutputStations):
  def __init__(self, **kwargs):
    _OutputStations.__init__(self, **kwargs)

  def make_plot(self, station, **kwargs):
    return _VelocityStations._make_plot(self, station, **kwargs)

  @staticmethod
  def from_netcdf(path):
    return _VelocityStations._from_netcdf(path)
