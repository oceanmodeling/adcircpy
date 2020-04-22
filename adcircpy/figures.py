import numpy as np
from matplotlib import rcParams
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


def get_topobathy_kwargs(values, vmin, vmax, colors=256):
    vmin = np.min(values) if vmin is None else vmin
    vmax = np.max(values) if vmax is None else vmax
    if vmax <= 0.:
        cmap = plt.cm.seismic
        col_val = 0.
        levels = np.linspace(vmin, vmax, colors)
    else:
        wet_count = int(np.floor(
            colors*(float((values < 0.).sum()) / float(values.size))))
        col_val = float(wet_count)/colors
        dry_count = colors - wet_count
        colors_undersea = plt.cm.bwr(np.linspace(1., 0., wet_count))
        colors_land = plt.cm.terrain(np.linspace(0.25, 1., dry_count))
        colors = np.vstack((colors_undersea, colors_land))
        cmap = LinearSegmentedColormap.from_list('cut_terrain', colors)
        wlevels = np.linspace(vmin, 0.0, wet_count, endpoint=False)
        dlevels = np.linspace(0.0, vmax, dry_count)
        levels = np.hstack((wlevels, dlevels))
    if vmax > 0:
        norm = FixPointNormalize(
            sealevel=0.0,
            vmax=vmax,
            vmin=vmin,
            col_val=col_val
            )
    else:
        norm = None
    return {'cmap': cmap,
            'norm': norm,
            'levels': levels,
            'col_val': col_val,
            # 'extend': 'both'
            }


def get_axes(axes, figsize=None, subplot=111):
    figsize = rcParams["figure.figsize"] if figsize is None else figsize
    if axes is None:
        fig = plt.figure(figsize=figsize)
        axes = fig.add_subplot(subplot)
    return axes


class FixPointNormalize(Normalize):
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


def _figure(f):
    def decorator(*argv, **kwargs):
        axes = get_axes(
            kwargs.get('axes', None),
            kwargs.get('figsize', None)
            )
        kwargs.update({'axes': axes})
        axes = f(*argv, **kwargs)
        if kwargs.get('show', False):
            axes.axis('scaled')
            plt.show()
        return axes
    return decorator
