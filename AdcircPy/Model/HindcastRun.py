from AdcircPy.Model._AdcircRun import _AdcircRun
from AdcircPy.Tides.TidalForcing import TidalForcing
from AdcircPy.Winds.BestTrack import BestTrack

class HindcastRun(_AdcircRun):
  def __init__(self, AdcircMesh, hurdat2_id, start_date=None, end_date=None, spinup_date=None, tides=True, netcdf=True, **kwargs):
    self.hurdat2_id   = hurdat2_id
    self.start_date   = start_date
    self.end_date     = end_date
    self.spinup_date  = spinup_date
    self.TidalForcing = tides
    self._init_BestTrack()
    self._init_TidalForcing()
    super(HindcastRun, self).__init__(AdcircMesh, **kwargs)
  
  def _init_BestTrack(self):
    self.BestTrack = BestTrack(self.hurdat2_id, self.start_date, self.end_date)

  def _init_TidalForcing(self):
    if self.TidalForcing == True:
      self.TidalForcing = TidalForcing(self.BestTrack.start_date, self.BestTrack.end_date, self.spinup_date)
    else:
      self.TidalForcing = None

  def _write_NWS(self):
    raise

  def _write_NRAMP(self):
    raise

  def _write_RNDAY(self):
    raise
  
  def _write_DRAMP(self):
    raise