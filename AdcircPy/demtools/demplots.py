from __future__ import absolute_import, division, print_function
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from AdcircPy.core._FixPointNormalize import FixPointNormalize

def plot_bathy(self, axes=None, vmin=None, vmax=None, title=None, colors=256):

    if axes is None:
        fig = plt.figure()
        axes = fig.add_subplot(111)

    if vmin is None:
        vmin = np.ceil(np.min(self.values))
    if vmax is None:
        vmax = np.floor(np.max(self.values))
    
    min_z = vmin
    max_z = vmax        

    if max_z < 0.:
        mlevel = np.mean(self.values)
        cmap = plt.cm.seismic
        col_val = 0.5
        levels=np.linspace(min_z, max_z, colors)        
    else:
        wet_count = int(np.floor(float(colors) * (float((self.values < 0.).sum())/float(self.values.size))))
        col_val = float(wet_count)/float(colors)
        dry_count = colors - wet_count
        colors_undersea = plt.cm.bwr(np.linspace(1., 0., wet_count))
        colors_land = plt.cm.terrain(np.linspace(0.25, 1., dry_count))
        colors = np.vstack((colors_undersea, colors_land))
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list('cut_terrain', colors)
        wlevels = np.linspace(min_z, 0.0, wet_count, endpoint=False)
        dlevels = np.linspace(0.0, max_z, dry_count)
        levels = np.hstack((wlevels, dlevels))
    norm = FixPointNormalize(sealevel=0.0, vmax=max_z, vmin=min_z, col_val=col_val)
    ax = axes.pcolormesh(self.x, np.flipud(self.y), self.values, cmap=cmap, norm=norm)
    axes.axis('scaled')
    if title is not None:
        axes.set_title(title)
    mappable = matplotlib.cm.ScalarMappable(cmap=cmap)
    mappable.set_array([])
    mappable.set_clim(min_z, max_z)
    divider = make_axes_locatable(axes)
    cax = divider.append_axes("bottom", size="2%", pad=0.5)
    cbar = plt.colorbar(mappable, cax=cax,
                        label=r'elevation [m]',
                        extend='both', orientation='horizontal')

    cbar.set_ticks([vmin, vmin + col_val *(vmax-vmin), vmax])
    cbar.set_ticklabels([np.around(vmin, 2), 0.0, np.around(vmax, 2)])
    
    return axes

def plot_diff(self, axes=None, vmin=None, vmax=None, title=None):
    if axes is None:
        fig = plt.figure()
        axes = fig.add_subplot(111)
    
    if vmin is None:
        vmin = np.ceil(np.min(self.values))
    if vmax is None:
        vmax = np.floor(np.max(self.values))

    cmap='seismic'
    norm = FixPointNormalize(sealevel=0, vmax=vmax, vmin=vmin, col_val=0.5)

    axes.imshow(self.values, origin='upper',
                                            extent=self.get_extent(),
                                            cmap=cmap,
                                            norm=norm)

    if title is not None:
        axes.set_title(title)

    mappable = matplotlib.cm.ScalarMappable(cmap=cmap)
    mappable.set_array([])
    mappable.set_clim(vmin, vmax)
    divider = make_axes_locatable(axes)
    cax = divider.append_axes("bottom", size="2%", pad=0.5)
    cbar = plt.colorbar(mappable, cax=cax,
                        label=r'elevation [$\Delta$m]',
                        extend='both', orientation='horizontal')

    cbar.set_ticks([vmin, vmin + 0.5*(vmax-vmin), vmax])
    cbar.set_ticklabels([np.around(vmin, 2), 0.0, np.around(vmax, 2)])
    return axes
    