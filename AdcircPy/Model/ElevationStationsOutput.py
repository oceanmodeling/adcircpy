from datetime import timedelta
from AdcircPy.Model._StationsOutput import _StationsOutput

class ElevationStationsOutput(_StationsOutput):
  def __init__(self, stations=dict(), sampling_frequency=timedelta(minutes=6), netcdf=True, spinup=False, harmonic_analysis=False):
    super(ElevationStationsOutput, self).__init__(sampling_frequency, netcdf, spinup, harmonic_analysis, **stations)
    self._comment1 = '! NOUTE,TOUTSE,TOUTFE,NSPOOLE:ELEV STATION OUTPUT INFO (UNIT  61)'
    self._comment2 = '! TOTAL NUMBER OF ELEVATION RECORDING STATIONS'
    
  @classmethod
  def from_csv(cls, path):
    raise NotImplementedError

  @classmethod
  def from_fort15(cls, path, sampling_frequency=None, netcdf=True, spinup=False, _hint='NOUTE', harmonic_analysis=False):
    return cls(cls._parse_fort15(path, _hint), sampling_frequency, netcdf, spinup, harmonic_analysis)


