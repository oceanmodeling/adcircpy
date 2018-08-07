import os
import numpy as np
import fnmatch
from netCDF4 import Dataset
from AdcircPy import Outputs
from AdcircPy.Outputs import _netcdf
from AdcircPy.Outputs import _ascii

def read_outputs(path, **kwargs):
    return Outputs.Outputs(path, **kwargs)._open_file()

def _open_file(self):
    if os.path.isfile(self._path)==False:
        raise FileNotFoundError("No such file or directory: %s" % path)
    self._check_netcdf()
    if self._nc == True:
        return self._read_netcdf()
    else:
        self._read_ascii_type()
        if self._ascii_type == 'harmonic_constituents':
            return Outputs.HarmonicConstituents.from_ascii(self._path, fort14=self._fort14, fort15=self._fort15, datum=self.datum, epsg=self.epsg)

# def _set__type(self):
#     if self._nc == True:
#         _netcdf._set__type(self)
#     else:
#         _ascii._set__type(self)