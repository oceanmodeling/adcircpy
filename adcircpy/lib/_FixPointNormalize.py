import numpy as np
from matplotlib.colors import Normalize


class _FixPointNormalize(Normalize):
    """
    This class is used for plotting. The reason it is declared here is that
    it is used by more than one submodule. In the future, this class will be
    native part of matplotlib. This definiton will be removed once the native
    matplotlib definition becomes available.

    Inspired by https://stackoverflow.com/questions/20144529/shifted-colorbar-matplotlib
    Subclassing Normalize to obtain a colormap with a fixpoint
    somewhere in the middle of the colormap.

    This may be useful for a `terrain` map, to set the "sea level"
    to a color in the blue/turquise range.
    """
    def __init__(self, vmin=None, vmax=None, sealevel=0, col_val=0.5,
                 clip=False):
        # sealevel is the fix point of the colormap (in data units)
        self.sealevel = sealevel
        # col_val is the color value in the range [0,1] that should represent
        # the sealevel.
        self.col_val = col_val
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.sealevel, self.vmax], [0, self.col_val, 1]
        if np.ma.is_masked(value)is False:
            value = np.ma.masked_invalid(value)
        return np.ma.masked_where(value.mask, np.interp(value, x, y))
