from AdcircPy.Tides.TidalForcing import TidalForcing
from AdcircPy import core

class TPXO(object):
  """  Runs at top level module import so that the TPXO cache gets initialized on disk. """
  core.init_TPXO_cache(core.get_cache_dir())
  @staticmethod
  def rebuild_cache():
    os.remove(core.get_cache_dir()+"/h_tpxo9.v1.nc")
    core.init_TPXO_cache(core.get_cache_dir())
