from AdcircPy.Model._AdcircRun import _AdcircRun
from AdcircPy.Tides import TidalForcing

class TidalRun(_AdcircRun):
  def __init__(self, AdcircMesh, start_date, end_date, spinup_date=None, constituents=None, **kwargs):
    self.TidalForcing = TidalForcing(start_date, end_date, spinup_date, constituents)
    super(TidalRun, self).__init__(AdcircMesh, **kwargs)
  
  def write_NWS(self):
    self.NWS=0
    self.f.write('{:<32d}'.format(self.NWS))

  def write_NRAMP(self):
    self.NRAMP=1
    self.f.write('{:<32d}'.format(self.NRAMP))

    