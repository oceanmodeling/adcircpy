from collections import OrderedDict
import os
import wget
import tarfile
from AdcircPy.Tides import _TidalDB

class TidalDB(OrderedDict):
  def __init__(self):
    self._init_constituents()
    self._units='rad/sec'


  def _init_constituents(self):
    _TidalDB._init_constituents(self)

  @staticmethod
  def init_TPXO_cache():
    _cachedir = os.getenv('LOCALAPPDATA')
    if _cachedir is None:
        _cachedir = os.getenv('HOME')+'/.cache/AdcircPy'
    else: 
        _cachedir += '/AdcircPy'
    os.makedirs(_cachedir, exist_ok=True)
    if os.path.isfile(_cachedir+"/h_tpxo9.v1.nc")==False:
        print('Building TPXO database cache on {}, please wait...'.format(_cachedir+"/h_tpxo9.v1.nc"))
        print('(This will only happen the first time you run this software)')
        url='ftp://ftp.oce.orst.edu/dist/tides/Global/tpxo9_netcdf.tar.gz'
        if os.path.isfile(_cachedir+"/tpxo9_netcdf.tar.gz")==False:
            tpxo=wget.download(url, out=_cachedir+"/tpxo9_netcdf.tar.gz")
            tpxo=tarfile.open(tpxo)
        else:
            tpxo=tarfile.open(tpxo)
        print(tpxo)