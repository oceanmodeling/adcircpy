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
  def from_fort15(cls, path, sampling_frequency=None, netcdf=True, spinup=False, harmonic_analysis=False, _hint='NOUTM'):
    return cls(cls.parse_fort15(path, _hint), sampling_frequency, netcdf, spinup, harmonic_analysis)