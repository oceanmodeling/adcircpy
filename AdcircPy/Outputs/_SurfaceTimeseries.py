import numpy as np
from AdcircPy.Surface import _fig

def make_animation(self, **kwargs)
    axes, idx = _fig.init_fig(self, axes, extent, title)
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
        
        