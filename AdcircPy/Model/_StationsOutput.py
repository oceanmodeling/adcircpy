from collections import OrderedDict
from AdcircPy.Model._ModelOutputs import _ModelOutputs

class _StationsOutput(OrderedDict, _ModelOutputs):
  def __init__(self, sampling_frequency, netcdf, spinup, harmonic_analysis, **kwargs):
    OrderedDict.__init__(self, **kwargs)
    _ModelOutputs.__init__(self, sampling_frequency, netcdf, spinup, harmonic_analysis)

  def add_station(self, station_id, x, y):
    self[station_id] = {'x' : x, 'y' : y}

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



