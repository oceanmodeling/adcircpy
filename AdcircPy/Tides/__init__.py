from collections import OrderedDict
from AdcircPy.Tides import TidalForcing as _TidalForcing


class TidalForcing(OrderedDict):
  """
    Notes:
    TPXO initializations is not part of TidalForcing because
    it depends on the mesh boundaries.
  """
  def __init__(self, start_date, end_date, constituents, spinup_date=None):
    self.constituents = constituents
    self.start_date   = start_date
    self.end_date     = end_date
    self.cachedir     = self.get_cache_dir()
    self.tpxo_path    = self.cachedir + "/h_tpxo9.v1.nc"
    self.init_constituent_dictionary()
    self.init_spinup_date(spinup_date)
    self.init_orbital_params()
    self.init_orbital_functions()
  
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

  def init_orbital_functions(self):
    _TidalForcing.init_orbital_functions(self)

  def _get_orbital_functions_start_of_record(self):
    return _TidalForcing._get_orbital_functions_start_of_record(self)

  def _get_orbital_functions_middle_of_record(self):
    return _TidalForcing._get_orbital_functions_middle_of_record(self)

  def _get_lunar_node_degrees(self, hours):
    return _TidalForcing._get_lunar_node_degrees(self, hours)

class TPXO(object):
  """  Runs at top level module import so that the TPXO cache gets initialized on disk. """
  _TidalForcing.init_TPXO_cache(_TidalForcing.get_cache_dir())
  @staticmethod
  def rebuild_cache():
    os.remove(_TidalForcing.get_cache_dir()+"/h_tpxo9.v1.nc")
    _TidalForcing.init_TPXO_cache(_TidalForcing.get_cache_dir())
