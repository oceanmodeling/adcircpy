from collections import OrderedDict
from AdcircPy.Tides import TidalForcing as _TidalForcing


class TidalForcing(OrderedDict):
  """
    Notes:
    TPXO initializations is not part of TidalForcing because
    it depends on the mesh boundaries.
  """
  def __init__(self, start_date, end_date, constituents=None, spinup_date=None):
    self.constituents = constituents
    self.start_date   = start_date
    self.end_date     = end_date
    self.cachedir     = self.get_cache_dir()
    self.tpxo_path    = self.cachedir + "/h_tpxo9.v1.nc"
    self.init_constituent_dictionary()
    self.init_spinup_date(spinup_date)
    self.init_orbital_params()
    self.init_node_factors()

  @property
  def units(self):
    return 'rad/sec'

  def init_constituent_dictionary(self):
    _TidalForcing.init_constituent_dictionary(self)

  @staticmethod
  def get_cache_dir():
    return _TidalForcing.get_cache_dir()

  def init_spinup_date(self, spinup_date):
    _TidalForcing.init_spinup_date(self, spinup_date)

  def init_orbital_params(self):
    _TidalForcing.init_orbital_params(self)

  def init_node_factors(self):
    _TidalForcing.init_node_factors(self)

  def _get_lunar_node(self, hours):
    return _TidalForcing._get_lunar_node(self, hours)

  def _get_lunar_perigee(self, hours):
    return _TidalForcing._get_lunar_perigee(self, hours)

  def _get_lunar_mean_longitude(self, hours):
    return _TidalForcing._get_lunar_mean_longitude(self, hours)

  def _get_solar_perigee(self, hours):
    return _TidalForcing._get_solar_perigee(self, hours)

  def _get_solar_mean_longitude(self, hours):
    return _TidalForcing._get_solar_mean_longitude(self, hours)

  def _get_nodal_factor(self, constituent):
    return _TidalForcing._get_nodal_factor(self, constituent)

  def _get_greenwich_term(self, constituent):
    return _TidalForcing._get_greenwich_term(self, constituent)

  def _EQ73(self):
    return _TidalForcing._EQ73(self)

  def _EQ74(self):
    return _TidalForcing._EQ74(self)

  def _EQ75(self):
    return _TidalForcing._EQ75(self)

  def _EQ76(self):
    return _TidalForcing._EQ76(self)

  def _EQ77(self):
    return _TidalForcing._EQ77(self)

  def _EQ78(self):
    return _TidalForcing._EQ78(self)

  def _EQ149(self):
    return _TidalForcing._EQ149(self)

  def _EQ197(self):
    return _TidalForcing._EQ197(self)

  def _EQ207(self):
    return _TidalForcing._EQ207(self)

  def _EQ213(self):
    return _TidalForcing._EQ213(self)
    
  def _EQ215(self):
    return _TidalForcing._EQ215(self)

  def _EQ227(self):
    return _TidalForcing._EQ227(self)

  def _EQ235(self):
    return _TidalForcing._EQ235(self)


class TPXO(object):
  """  Runs at top level module import so that the TPXO cache gets initialized on disk. """
  _TidalForcing.init_TPXO_cache(_TidalForcing.get_cache_dir())
  @staticmethod
  def rebuild_cache():
    os.remove(_TidalForcing.get_cache_dir()+"/h_tpxo9.v1.nc")
    _TidalForcing.init_TPXO_cache(_TidalForcing.get_cache_dir())
