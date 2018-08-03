import os
import numpy as np
import fnmatch
from netCDF4 import Dataset
from AdcircPy import Outputs
from AdcircPy.Outputs import _netcdf
from AdcircPy.Outputs import _ascii

def read_outputs(path, **kwargs):
    return Outputs.Outputs(path)._open_file()

def _open_file(self):
    if os.path.isfile(self._path)==False:
        raise FileNotFoundError("No such file or directory: %s" % path)
    self._check_netcdf()
    if self._nc == True:
        return self._read_netcdf()
    else:
        return self._read_ascii()

# def _set__type(self):
#     if self._nc == True:
#         _netcdf._set__type(self)
#     else:
#         _ascii._set__type(self)