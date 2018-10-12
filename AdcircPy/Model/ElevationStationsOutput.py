from AdcircPy.Model._AdcircOutputs import _AdcircOutputs

class ElevationStationsOutput(_AdcircOutputs):
  def __init__(self, stations, sampling_frequency, netcdf, spinup):
    super(ElevationStationsOutput, self).__init__(sampling_frequency, netcdf, spinup, **stations)

  @classmethod
  def from_csv(cls, path):
    raise NotImplementedError

  @classmethod
  def from_fort15(cls, path, sampling_frequency=None, netcdf=True, spinup=False):
    if path is None:
      return
    stations=dict()
    with open(path,'r') as f:
      for line in f:
        if 'NOUTE' in line:
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
    return cls(stations, sampling_frequency, netcdf, spinup)

