import matplotlib.pyplot as plt
import numpy as np
from AdcircPy import fig
from AdcircPy import Mesh

def get_difference(self, other):
    params = {  'x'                   : self.x,
                'y'                   : self.y,
                'elements'            : self.elements,
                'values'              : self.values - other.values,
                'nodeID'              : self.nodeID,
                'elementID'           : self.elementID, 
                "ocean_boundaries"    : self.ocean_boundaries,
                "land_boundaries"     : self.land_boundaries,
                "inner_boundaries"    : self.inner_boundaries,
                "weir_boundaries"     : self.weir_boundaries,
                "inflow_boundaries"   : self.inflow_boundaries,
                "outflow_boundaries"  : self.outflow_boundaries,
                "culvert_boundaries"  : self.culvert_boundaries}
    return Mesh.SurfaceDifference(**params)

def plot_diff(self, extent=None, axes=None, vmin=None, vmax=None, title=None, **kwargs):
    axes, idx = fig._init_fig(self, axes, extent, title)
    if vmin is None:
        vmin = np.min(self.values[idx])
    if vmax is None:
        np.max(self.values[idx])

    cmap = plt.get_cmap(kwargs.pop("cmap", "seismic"))
    levels = kwargs.pop("levels", np.linspace(np.min(self.values[idx]), np.max(self.values[idx]), 256))
    norm = fig.FixPointNormalize(sealevel=0, vmax=vmax, vmin=vmin, col_val=0.5)
    if np.ma.is_masked(self.values):
        trimask = np.any(self.values.mask[self.elements], axis=1)
        Tri = matplotlib.tri.Triangulation(self.x, self.y, self.elements, trimask)
        axes.tricontourf(Tri, self.values, levels=levels, cmap=cmap, extend='both', norm=norm)
    else:
        axes.tricontourf(self.x, self.y, self.elements, self.values, levels=levels, cmap=cmap, extend='both', norm=norm)
    cbar = fig._init_colorbar(axes, cmap, vmin, vmax)
    cbar.set_ticks([vmin, vmin + 0.5*(vmax-vmin), vmax])
    cbar.set_ticklabels([np.around(vmin, 2), 0.0, np.around(vmax, 2)])
    cbar.set_label(r'elevation [$\Delta$ m]')
    return axes