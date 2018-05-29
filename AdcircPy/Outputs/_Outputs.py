import os
from netCDF4 import Dataset
from AdcircPy.Outputs import _netcdf
from AdcircPy.Outputs import _ascii

def read_outputs(path, **kwargs):
    # Keep the returns outside of the try/except in case of throwbacks.
    if os.path.isfile(path)==False:
        raise FileNotFoundError("No such file or directory: %s" % path)
    try:
        nc = Dataset(path)
        nc = True
    except:
        nc = False
    if nc == True:
        return _netcdf.read_netcdf_output(path, **kwargs)
    elif nc == False:
        return _ascii.read_ascii_output(path, **kwargs)