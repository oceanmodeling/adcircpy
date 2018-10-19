from AdcircPy.Model._StationsOutput import _StationsOutput

class _HarmonicConstituentsStations(dict):
  def __init__(self, fort15, _hint):
    super(_HarmonicConstituentsStations, self).__init__()
    self._fort15 = fort15
    self._hint = _hint
    self._init_fort15()

  def _init_fort15(self):
    if self._fort15 is None:
      raise Exception('The fort.15 file used to generate this output is required for post processing.')
    else:
      self._fort15 = _StationsOutput.parse_fort15(self._fort15, self._hint)

