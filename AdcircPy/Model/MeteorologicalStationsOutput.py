from datetime import timedelta
from AdcircPy.Model._StationsOutput import _StationsOutput

class MeteorologicalStationsOutput(_StationsOutput):
  def __init__(self, stations=dict(), sampling_frequency=timedelta(minutes=6), netcdf=True, spinup=False, harmonic_analysis=False):
    super(MeteorologicalStationsOutput, self).__init__(sampling_frequency, netcdf, spinup, harmonic_analysis, **stations)
    self._comment1 = '! NOUTM,TOUTSM,TOUTFM,NSPOOLM:VEL STATION OUTPUT INFO (UNIT  62)'
    self._comment2 = '! TOTAL NUMBER OF VELOCITY RECORDING STATIONS'
  
  @classmethod
  def from_csv(cls, path):
    raise NotImplementedError

  @classmethod
  def from_fort15(cls, path, sampling_frequency=None, netcdf=True, spinup=False, harmonic_analysis=False):
    stations=dict()
    with open(path,'r') as f:
      for line in f:
        if 'NOUTM' in line:
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
    return cls(stations, sampling_frequency, netcdf, spinup, harmonic_analysis)