from AdcircPy.Model._AdcircRun import _AdcircRun
from AdcircPy.Tides import TidalForcing

class TidalRun(_AdcircRun):
  def __init__(self, AdcircMesh, start_date, end_date,
               spinup_date=None, constituents=None, **kwargs):
    self.TidalForcing = TidalForcing(start_date, end_date, spinup_date, constituents)
    super(TidalRun, self).__init__(AdcircMesh, **kwargs)
 
  def _init_NWS(self):
    self.NWS=0

  def _init_DRAMP(self):
    if self.DRAMP is None:
      self.DRAMP = ((2/3)*(self.TidalForcing.start_date - self.TidalForcing.spinup_date).total_seconds())/(60*60*24)

  def _write_NRAMP(self):
    self.NRAMP=1
    self.f.write('{:<32d}'.format(self.NRAMP))

  def _write_RNDAY(self):
    if self.IHOT==0:
      RNDAY = (self.TidalForcing.start_date - self.TidalForcing.spinup_date).total_seconds()/(60*60*24)
    elif self.IHOT==567:
      RNDAY = (self.TidalForcing.end_date - self.TidalForcing.spinup_date).total_seconds()/(60*60*24)
    self.f.write('{:<32.2f}'.format(RNDAY))
  
  def _write_DRAMP(self):
    if self.DRAMP is None:
      self.DRAMP = ((2/3)*(self.TidalForcing.start_date - self.TidalForcing.spinup_date).total_seconds())/(60*60*24)
    self.f.write('{:<32.1f}'.format(self.DRAMP))


    
