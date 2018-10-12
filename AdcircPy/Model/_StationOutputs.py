from collections import OrderedDict
from datetime import timedelta

class _StationOutputs(OrderedDict):
  def __init__(self, sampling_frequency, netcdf, spinup, **kwargs):
    super(_StationOutputs, self).__init__(**kwargs)
    self._init_sampling_frequency(sampling_frequency)
    self.spinup = spinup
    self.netcdf = netcdf

  def add_station(self, station_id, x, y):
    self[station_id] = {'x' : x, 'y' : y}

  def _init_sampling_frequency(self, sampling_frequency):
    if sampling_frequency is None:
      self.sampling_frequency=timedelta(minutes=6)
    elif isinstance(sampling_frequency, timedelta):
      self.sampling_frequency=sampling_frequency
    else:
      raise TypeError('sampling_frequency must be a datetime.timedelta object.')

  @staticmethod
  def _parse_fort15(path, _hint):
    stations=dict()
    with open(path,'r') as f:
      for line in f:
        if _hint in line:
          num = int(f.readline().split('!')[0].strip().split(' ')[0])
          for i in range(num):
            line = f.readline().split('!')
            if len(line)>0:
              station_name = line[1].strip(' \n')
            else:
              station_name = str(i)
            coords = line[0].split(' ')
            coords = [float(x) for x in coords if x != '']
            x = float(coords[0])
            y = float(coords[1])
            stations[station_name]={'x' : x, 'y' : y}
    return stations



