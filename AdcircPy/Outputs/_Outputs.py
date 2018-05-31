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
        
        
def surface_animation(self, extent=None, axes=None, title=None, **kwargs):
    axes, idx = fig._init_fig(self, axes, extent, title)
    start_slice = kwargs.pop("start_slice", 0)
    stop_slice  = kwargs.pop("stop_slice", len(self.timestep))
    vals = list()
    for values in self.values[start_slice:stop_slice]:
        vals.append(values[idx])
    vals = np.asarray(vals)
    vals = np.ma.masked_where(vals==-99999.0, vals)
    vmin = kwargs.pop("vmin", np.min(vals))
    vmax = kwargs.pop("vmax", np.max(vals))
    return axes