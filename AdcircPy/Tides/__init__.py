from AdcircPy.Tides.TidalForcing import TidalForcing
from AdcircPy.Tides import _TPXO

class TPXO(object):
  """  Runs at top level module import so that the TPXO cache gets initialized on disk. """
  _TPXO.init_TPXO_cache(_TPXO.get_cache_dir())
  @staticmethod
  def rebuild_cache():
    os.remove(_TPXO.get_cache_dir()+"/h_tpxo9.v1.nc")
    _TPXO.init_TPXO_cache(_TPXO.get_cache_dir())
