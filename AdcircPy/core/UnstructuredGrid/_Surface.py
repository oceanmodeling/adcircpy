from collections import defaultdict
import numpy as np
from matplotlib.tri import Triangulation
import matplotlib.pyplot as plt
from matplotlib.path import Path
from scipy.interpolate import griddata
import pyproj
from AdcircPy import Surface
from AdcircPy import demtools
from AdcircPy.Surface import _fig

def make_plot(self, extent=None, epsg=None, axes=None, title=None, step=0, total_colors=256, cbar_label=None, **kwargs):
    axes, idx = _fig.init_fig(self, axes, extent, title, epsg)
    if isinstance(self.values, list):
        values = self.values[step]
    else:
        values = self.values
    vmin = kwargs.pop("vmin", np.min(values))
    vmax = kwargs.pop("vmax", np.max(values))
    cmap = kwargs.pop("cmap", "jet")
    levels = kwargs.pop("levels", np.linspace(vmin, vmax, total_colors))
    if np.ma.is_masked(values):
        trimask = np.any(values.mask[self.elements], axis=1)
        Tri = Triangulation(self.x, self.y, self.elements, trimask)
        axes.tricontourf(Tri, values, levels=levels, cmap=cmap, extend='both')
    else:
        axes.tricontourf(self.x, self.y, self.elements, values, levels=levels, cmap=cmap, extend='both')
    cbar = _fig.init_colorbar(axes, cmap, vmin, vmax)
    if cbar_label is not None:
        cbar.set_label(cbar_label)
    cbar.set_ticks([vmin,
                    vmin+(1./4.)*(vmax-vmin),
                    vmin+(1./2.)*(vmax-vmin),
                    vmin+(3./4.)*(vmax-vmin),
                    vmax])
    cbar.set_ticklabels([np.around(vmin, 2), 
                            np.around(vmin+(1./4.)*(vmax-vmin), 2),
                            np.around(vmin+(1./2.)*(vmax-vmin), 2),
                            np.around(vmin+(3./4.)*(vmax-vmin), 2),
                            np.around(vmax, 2)])
    return axes

def plot_velocity(self, **kwargs):
    raise NotImplementedError("Coming soon!")
    start_timestep = kwargs.pop('start_timestep', 0)
    stop_timestep  = kwargs.pop('stop_timestep', len(self.time))
    if axes is None:                
        fig = plt.figure()
        axes  = fig.add_subplot(111)
    ax = axes.quiver(self.x, self.y, self.u[0,:], self.v[:,0], vmin=vmin, vmax=vmax, **kwargs)
    axes.axis('scaled')
    def update(i):
       ax.set_array(self.Dataset['zs'][i,:-1,:-1].ravel())
       return ax
    
    anim = FuncAnimation(fig, update, frames=np.arange(start_timestep+1, stop_timestep), interval=interval)
    if colorbar==True:
        plt.colorbar(ax)
    if show is True:
        plt.show()
    return anim

def get_dict(self):
    return {  'x'                   : self.x,
              'y'                   : self.y,
              'elements'            : self.elements,
              'values'              : self.values,
              'nodeID'              : self.nodeID,
              'elementID'           : self.elementID, 
              "ocean_boundaries"    : self.ocean_boundaries,
              "land_boundaries"     : self.land_boundaries,
              "inner_boundaries"    : self.inner_boundaries,
              "weir_boundaries"     : self.weir_boundaries,
              "inflow_boundaries"   : self.inflow_boundaries,
              "outflow_boundaries"  : self.outflow_boundaries,
              "culvert_boundaries"  : self.culvert_boundaries}


def get_difference(self, other, step=0):
    

    if isinstance(self.values, list):
        self_values = self.values[step]
    else:
        self_values = self.values

    if isinstance(other.values, list):
        other_values = other.values[step]
    else:
        other_values = other.values

    if np.ma.is_masked(self_values):
        self_values = np.ma.filled(self_values, 0.)

    if np.ma.is_masked(other.values):
        other_values = np.ma.filled(other_values, 0.)
    values = self_values - other_values
    # values = np.ma.masked_equal(values, 0.)
    kwargs = self.get_dict()
    kwargs['values'] = values
    return Surface.SurfaceDifference(**kwargs)

def get_mean_value(self, **kwargs):
    epsg   = kwargs.pop("epsg", self.epsg)
    extent = kwargs.pop("extent", self.get_extent(epsg=epsg))
    step   = kwargs.pop("step", 0)
    idx = self.get_extent_idx(extent, epsg)
    if isinstance(self.values, list):
        values = self.values[step]
    else:
        values = self.values
    return np.nanmean(values[idx])

def plot_diff(self, extent=None, epsg=None, axes=None, vmin=None, vmax=None, title=None, **kwargs):
    axes, idx = _fig.init_fig(self, axes, extent, title, epsg)
    if vmin is None:
        vmin = np.min(self.values[idx])
    if vmax is None:
        np.max(self.values[idx])

    cmap = plt.get_cmap(kwargs.pop("cmap", "seismic"))
    levels = kwargs.pop("levels", np.linspace(np.min(self.values[idx]), np.max(self.values[idx]), 256))
    norm = _fig.FixPointNormalize(sealevel=0, vmax=vmax, vmin=vmin, col_val=0.5)
    if np.ma.is_masked(self.values):
        trimask = np.any(self.values.mask[self.elements], axis=1)
        Tri = Triangulation(self.x, self.y, self.elements, trimask)
        axes.tricontourf(Tri, self.values, levels=levels, cmap=cmap, extend='both', norm=norm)
    else:
        axes.tricontourf(self.x, self.y, self.elements, self.values, levels=levels, cmap=cmap, extend='both', norm=norm)
    cbar = _fig.init_colorbar(axes, cmap, vmin, vmax)
    cbar.set_ticks([vmin, vmin + 0.5*(vmax-vmin), vmax])
    cbar.set_ticklabels([np.around(vmin, 2), 0.0, np.around(vmax, 2)])
    cbar.set_label(r'elevation [$\Delta$ m]')
    return axes

def rasterize_to_geoTransform(self, geoTransform, shape, **kwargs):
    """
    Converts an ADCIRC mesh object into a ADCIRC raster object.
    With this

    Parameters
    ----------
    geoTransform : tuple
        geoTransform[0] /* top left x */
        geoTransform[1] /* west-east pixel resolution */
        geoTransform[2] /* 0 */
        geoTransform[3] /* top left y */
        geoTransform[4] /* 0 */
        geoTransform[5] /* north-south pixel resolution (negative value) */

    xpixels : int
        Description of arg2

    ypixels : int

    Returns
    -------
    raster object
        Can be used to export to geoTif, and other operations such as filtering.
    """
    epsg = kwargs.pop("epsg", self.epsg)
    padding = kwargs.pop("padding", None)
    xpixels, ypixels = shape
    x = np.linspace(geoTransform[0], geoTransform[0] + xpixels*geoTransform[1], xpixels)
    y = np.linspace(geoTransform[3] + ypixels*geoTransform[5], geoTransform[3], ypixels)
    xt, yt = np.meshgrid(x, y)
    xt = xt.reshape(xt.size)
    yt = np.flipud(yt.reshape(yt.size))
    xyt = np.vstack((xt,yt)).T
    # create path object of target bounding box
    bbox_path = Path([(np.min(xt), np.min(yt)),
                    (np.max(xt), np.min(yt)),
                    (np.max(xt), np.max(yt)),
                    (np.min(xt), np.max(yt)),
                    (np.min(xt), np.min(yt))], closed=True)
    # interpolate mesh information to bounding box grid.
    xyz = self.get_xyz(extent=bbox_path)
    zt = griddata((xyz[:,0], xyz[:,1]), xyz[:,2], (xt, yt), method='linear', fill_value=np.nan)
    # Generate boundary masks.
    outerBoundary = self.build_outer_polygon()
    outerBoundary = outerBoundary.clip_to_bbox((np.min(xt), np.min(yt), np.max(xt), np.max(yt)))
    mask = np.logical_or(np.isin(zt, np.nan), ~outerBoundary.contains_points(xyt))
    zt = np.ma.masked_array(zt, mask)
    innerBoundaries = self.build_inner_polygons()
    for innerBoundary in innerBoundaries:
        if outerBoundary.intersects_path(innerBoundary):
            innerBoundary = innerBoundary.clip_to_bbox((np.min(xt), np.min(yt), np.max(xt), np.max(yt)))
            mask = np.logical_or(zt.mask, innerBoundary.contains_points(xyt))
            zt = np.ma.masked_array(zt.data, mask)
    #TODO: Compound mask for other boundaries.
    if self.weir_boundaries is not None:
        for boundary in self.weir_boundaries:
            pass
    if self.culvert_boundaries is not None:
        for boundary in self.culvert_boundaries:
            pass
    if padding is not None:
        padding = griddata((padding[:,0], padding[:,1]), padding[:,2], (xt, yt), method='nearest')
        idx, = np.where(zt.mask)
        zt = np.ma.filled(zt, np.nan)
        zt[idx] = padding[idx]
    zt = zt.reshape(shape)
    return demtools.DEM(x, y, zt, geoTransform, self.epsg, self.datum)

def get_raster_from_extent(self, extent , dx, dy, epsg, padding=None):
    """
    This script rasterizes the mesh into a regular grid that can be exported as GeoTiff.
    Uselful for exploring the data on a GIS program.
    """
    min_x, max_x, min_y, max_y = extent
    
    x = np.arange(min_x, max_x+dx, dy)
    y = np.arange(min_y, max_y+dy, dy)
    
    xt, yt = np.meshgrid(x, y)
    shape = xt.shape
    
    xt = xt.reshape(xt.size)
    yt = np.flipud(yt.reshape(yt.size))
    xyt = np.vstack((xt,yt)).T

    bbox_path = Path([(np.min(xt), np.min(yt)),
                    (np.max(xt), np.min(yt)),
                    (np.max(xt), np.max(yt)),
                    (np.min(xt), np.max(yt)),
                    (np.min(xt), np.min(yt))], closed=True)

    xyz = self.get_xyz(extent=bbox_path, radius=0.2)
    zt = griddata((xyz[:,0], xyz[:,1]), xyz[:,2], (xt, yt), method='linear', fill_value=np.nan)
    
    
    # generate masks
    if ~hasattr(self,'outer_boundary'):
        self.outer_boundary = self.build_outer_polygon()
    path = self.outer_boundary.clip_to_bbox((np.min(xt), np.min(yt), np.max(xt), np.max(yt)))
    mask1 = np.isin(zt, -99999.0)
    mask2 = path.contains_points(xyt)
    mask = np.logical_or(mask1, ~mask2)
    for i in range(len(self.land_boundaries)):
        if self.land_boundaries[i][-1] in [1,21]:
            extents = self.land_boundaries[i][0].get_extents()
            extents = extents.get_points()
            c1 = np.logical_and(np.min(xt) <= extents[1,0], np.max(xt) >= extents[0,0])
            c2 = np.logical_and(np.min(yt) <= extents[1,1], np.max(yt) >= extents[0,1])
            if np.logical_and(c1,c2):
                shape = self.land_boundaries[i][0].clip_to_bbox((np.min(xt), np.min(yt), np.max(xt), np.max(yt)))
                mask3 = shape.contains_points(xyt)
                mask = np.logical_or(mask, mask3)  
    geoTransform = (np.min(min_x), dx, 0, np.max(y), 0, -np.abs(dy))
    depth = np.ma.masked_array(zt, mask).reshape(_shape)
    if padding is not None:
        depth = adcpy.adcirc.raster._apply_padding(x, y, depth, padding)
    DEM = demtools.DEM()
    return DEM(x, y, depth, geoTransform, 4326, self.datum)


def get_contours(self, levels, **kwargs):
    epsg = kwargs.pop("epsg", self.epsg)

    if epsg != self.epsg:
        self_proj   = pyproj.Proj(init='epsg:{}'.format(self.epsg))
        target_proj = pyproj.Proj(init='epsg:{}'.format(epsg))
        x, y = pyproj.transform(self_proj, target_proj, self.x, self.y)

    else:
        x = self.x
        y = self.y

    if np.ma.is_masked(self.values):
        mask = np.where(self.values.mask)
        trimask = np.any(np.in1d(self.elements, mask).reshape(-1, 3), axis=1)
        Tri = Triangulation(x, y, self.elements, trimask)
    else:
        Tri = Triangulation(x, y, self.elements)
    fig = plt.figure()
    axes = fig.add_subplot(111)
    ax = axes.tricontour(Tri, self.values, levels=levels)
    plt.close(fig)

    contours = defaultdict(list)
    for i, LineCollection in enumerate(ax.collections):
        for Path in LineCollection.get_paths():
            contours[ax.levels[i]].append(Path)
    return dict(contours)

def get_values_at_xy(self, x, y, step=0, method='linear'):
    if isinstance(self.values, list):
        values = self.values[step]
    else:
        values = self.values
    if np.ma.is_masked(values):
        values = np.ma.filled(values, 0.0)
    elif np.isin(values, -99999.0).any():
        idx = np.where(np.isin(values, -99999.0))
        values[idx] = 0.0    
    if method != 'force':
        return griddata((self.x,self.y),values,(x,y), method=method)
    else:
        idx = np.where(~np.isin(values, 0.0))
        return griddata((self.x[idx], self.y[idx]), values[idx], (x, y), method='nearest')
