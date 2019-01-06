from AdcircPy.Model._AdcircRun import _AdcircRun
from AdcircPy.Tides.TidalForcing import TidalForcing
from AdcircPy.Winds.BestTrack import BestTrack

class BestTrackRun(_AdcircRun):
  def __init__(self, AdcircMesh, storm_id, start_date=None, end_date=None, spinup_date=None, tides=True, netcdf=True, **kwargs):
    self.storm_id   = storm_id
    self.start_date   = start_date
    self.end_date     = end_date
    self.spinup_date  = spinup_date
    self.TidalForcing = tides
    self._init_fort22()
    self._init_TidalForcing(kwargs.pop('constituents', None))
    super(BestTrackRun, self).__init__(AdcircMesh, **kwargs)
  
  def _init_fort22(self):
    self.fort22 = BestTrack(self.storm_id, self.start_date, self.end_date)

  def _init_TidalForcing(self, constituents):
    if self.TidalForcing == True:
      self.TidalForcing = TidalForcing(self.fort22.start_date, self.fort22.end_date, self.spinup_date, constituents=constituents)
    else:
      self.TidalForcing = None

  def _init_NWS(self):
    self.NWS=20


  def _init_DRAMP(self):
    if self.DRAMP is None:
      if self.TidalForcing is not None:
        self.DRAMP = ((2/3)*(self.TidalForcing.start_date - self.TidalForcing.spinup_date).total_seconds())/(60*60*24)
      else:
        self.DRAMP = 0


  def _write_NRAMP(self):
    if self.IHOT>0:
      self.NRAMP=8
      self.f.write('{:<32d}'.format(self.NRAMP))
    else:
      self.NRAMP=0
      self.f.write('{:<32d}'.format(self.NRAMP))
    
  def _write_RNDAY(self):
    if self.TidalForcing is not None:
      if self.IHOT==0:
        RNDAY = (self.TidalForcing.start_date - self.TidalForcing.spinup_date).total_seconds()/(60*60*24)
      elif self.IHOT==567:
        RNDAY = (self.TidalForcing.end_date - self.TidalForcing.spinup_date).total_seconds()/(60*60*24)
    else:
      RNDAY = (self.fort22.end_date - self.fort22.start_date).total_seconds()/(60*60*24)
    self.f.write('{:<32.2f}'.format(RNDAY))
  

  def _write_DRAMP(self):
    if self.NRAMP>0 and self.TidalForcing is not None:
      if self.DUnRampMete is None:
        self.DUnRampMete = (self.TidalForcing.start_date - self.TidalForcing.spinup_date).days
      else:
        self.DUnRampMete = 0
      if self.DRAMPElev is None:
        self.DRAMPElev = self.DRAMP
      if self.DRAMPTip is None:
        self.DRAMPTip = self.DRAMP
      self.f.write('{:<4.1f} '.format(self.DRAMP))
      self.f.write('{:<2}'.format(int(self.DRAMPExtFlux)))
      self.f.write('{:<2}'.format(int(self.FluxSettlingTime)))
      self.f.write('{:<2}'.format(int(self.DRAMPIntFlux)))
      self.f.write('{:<3}'.format(int(self.DRAMPElev)))
      self.f.write('{:<3}'.format(int(self.DRAMPTip)))
      self.f.write('{:<4.1f}'.format(self.DRAMPMete))
      self.f.write('{:<2}'.format(int(self.DRAMPWRad)))
      self.f.write('{:<2}'.format(int(self.DUnRampMete)))
      self.f.write('{:<7}'.format(''))
    elif self.NRAMP==0 and self.TidalForcing is not None:
      self.DRAMP = ((2/3)*(self.TidalForcing.start_date - self.TidalForcing.spinup_date).total_seconds())/(60*60*24)
      self.f.write('{:<32.1f}'.format(self.DRAMP))
    else:
      raise NotImplementedError('I suppose this is a met-only run. We need to code more.')

  def dump(self, path):
    """
    Overloads the parent's dump method such that if there is tidal forcing present
    it can call the parent, otherwise the met only run requires a slightly different
    fort.15 generation as met-only has no coldstart phase.
    """
    if self.TidalForcing is not None:
      _AdcircRun.dump(self, path)
    else:
      self.IHOT=0
      self.NRAMP=8
      self._write_fort15()
    self.fort22.dump(path)
