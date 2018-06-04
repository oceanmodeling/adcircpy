import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable

class FixPointNormalize(Normalize):
    """ 
    Inspired by https://stackoverflow.com/questions/20144529/shifted-colorbar-matplotlib
    Subclassing Normalize to obtain a colormap with a fixpoint 
    somewhere in the middle of the colormap.
    This may be useful for a `terrain` map, to set the "sea level" 
    to a color in the blue/turquise range.
    """
    def __init__(self, vmin=None, vmax=None, sealevel=0, col_val = 0.5, clip=False):
        # sealevel is the fix point of the colormap (in data units)
        self.sealevel = sealevel
        # col_val is the color value in the range [0,1] that should represent the sealevel.
        self.col_val = col_val
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.sealevel, self.vmax], [0, self.col_val, 1]
        return np.ma.masked_where(value.mask, np.interp(value, x, y))

def init_fig(self, axes, extent, title):
    if axes is None:                
        fig = plt.figure()
        axes  = fig.add_subplot(111)
    if title is not None:
        axes.set_title(title)
    if extent is None:
        extent = self.get_extent()
    idx = self.get_extent_idx(extent=extent)
    axes.axis('scaled')
    axes.axis(extent) 
    return axes, idx

def init_colorbar(axes, cmap, vmin, vmax):
    mappable = ScalarMappable(cmap=cmap)
    divider = make_axes_locatable(axes)
    cax = divider.append_axes("bottom", size="2%", pad=0.5)
    mappable.set_array([])
    mappable.set_clim(vmin, vmax)
    return plt.colorbar(mappable, cax=cax, extend='both', orientation='horizontal')