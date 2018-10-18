import os
import wget
import tarfile
from datetime import datetime, timedelta
from AdcircPy.core.FixPointNormalize import FixPointNormalize
"""
Function for aliasing kwargs used in _AdcircRun.py
Reference:
https://stackoverflow.com/questions/29374425/aliasing-in-the-names-of-function-arguments-double-naming
"""
import functools
def alias(aliases):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*argv, **kwargs):
            for name, alias in aliases.items():
                if name not in kwargs and alias in kwargs:
                    kwargs[name] = kwargs[alias]
            return func(*argv, **kwargs)
        return wrapper
    return decorator

"""
Functions for initializing the TPXO cache
"""

def init_TPXO_cache(cachedir):
  os.makedirs(cachedir, exist_ok=True)
  if os.path.isfile(cachedir+"/h_tpxo9.v1.nc")==False:
    print('Building TPXO database cache on {}, please wait...'.format(cachedir+"/h_tpxo9.v1.nc"))
    print('(This will only happen the first time you run this software)')
    url='ftp://ftp.oce.orst.edu/dist/tides/Global/tpxo9_netcdf.tar.gz'
    if os.path.isfile(cachedir+"/tpxo9_netcdf.tar.gz")==False:
      wget.download(url, out=cachedir+"/tpxo9_netcdf.tar.gz")
    tpxo=tarfile.open(cachedir+"/tpxo9_netcdf.tar.gz")
    tpxo.extract('h_tpxo9.v1.nc', path=cachedir) # untar file into same directory
    tpxo.close()
    os.remove(cachedir+"/tpxo9_netcdf.tar.gz")

def get_cache_dir():
  cachedir = os.getenv('LOCALAPPDATA')
  if cachedir is None:
    cachedir = os.getenv('HOME')+'/.cache/AdcircPy'
  else: 
    cachedir += '/AdcircPy'
  return cachedir
