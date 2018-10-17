
class _OutputStations(dict):
  def __init__(self, time, **station_data):
    self.time = time
    dict.__init__(self, **station_data)
