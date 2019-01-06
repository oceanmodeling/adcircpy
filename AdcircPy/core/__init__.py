import os
import wget
import tarfile
from datetime import datetime, timedelta
from AdcircPy.core._Figure import _Figure
from AdcircPy.core._FixPointNormalize import _FixPointNormalize
from AdcircPy.core._PostGIS import _PostGIS
from AdcircPy.core.PBS import PBS
from AdcircPy.core.ServerConfiguration import ServerConfiguration

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
