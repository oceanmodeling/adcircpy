import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import griddata
from haversine import haversine
from AdcircPy.core._FixPointNormalize import FixPointNormalize
from AdcircPy.Mesh import _fort14
from AdcircPy.Mesh import _fort13

def _from_fort14(cls, fort14, datum='MSL', epsg=4326, fort13=None, fort15=None, datum_grid=None):
    kwargs = _fort14.parse_fort14(fort14)
    kwargs['datum'] = datum
    kwargs['epsg']  = epsg
    kwargs['fort13'] = fort13
    kwargs['fort15'] = fort15
    kwargs['datum_grid'] = datum_grid
    return cls(**kwargs)

def _make_plot(self, extent=None, epsg=None, axes=None, title=None, total_colors=256, cbar_label=r'elevation [m]', **kwargs):
    self._init_fig(axes, extent, title, epsg)
    vmin = kwargs.pop("vmin", np.min(self.values[self._idx]))
    vmax = kwargs.pop("vmax", np.max(self.values[self._idx]))
    show = kwargs.pop("show", False)
    if vmax < 0.:

        levels = kwargs.pop("levels", np.linspace(vmin, vmax, total_colors))
        mlevel = np.mean(self.values[self._idx])
        cmap   = kwargs.pop("cmap", plt.cm.seismic)
        col_val = 0.5      
    else:
        mlevel = 0.
        wet_count = int(np.floor(total_colors * \
                (float((self.values[self._idx] < 0.).sum())/self.values[self._idx].size)))
        col_val = float(wet_count)/float(total_colors)
        dry_count = total_colors - wet_count
        colors_undersea = plt.cm.bwr(np.linspace(1., 0., wet_count))
        colors_land = plt.cm.terrain(np.linspace(0.25, 1., dry_count))
        colors = np.vstack((colors_undersea, colors_land))
        cmap = LinearSegmentedColormap.from_list('cut_terrain', colors)
        wlevels = np.linspace(vmin, mlevel, wet_count, endpoint=False)
        dlevels = np.linspace(mlevel, vmax, dry_count)
        levels = np.hstack((wlevels, dlevels))
    norm = FixPointNormalize(sealevel=mlevel, vmax=vmax, vmin=vmin, col_val=col_val)
    self._axes.tricontourf(self.x, self.y, self.elements, self.values,
                    levels=levels, cmap=cmap, norm=norm, extend='both', **kwargs)
    self._init_cbar(self._axes, cmap, vmin, vmax)
    self._cbar.set_label(cbar_label)
    self._cbar.set_ticks([vmin, vmin + col_val *(vmax-vmin), vmax])
    self._cbar.set_ticklabels([np.around(vmin, 2), mlevel, np.around(vmax, 2)])
    if show == True:
        plt.show()
    return self._axes

def plot_shoreline(self, extent=None, epsg=None, axes=None, title=None, color='black', linewidth=0.5):
    axes, idx = Figures._init_fig(self, axes, extent, title, epsg)
    if type(self.values) is list:
        if self.bathymetry is None:
            raise TypeError("fort.14 required to plot shoreline")
        values = self.bathymetry
    else:
        values = self.values

    axes.tricontour(self.x, self.y, self.elements, values, levels=[0.], colors=color, linewidths=linewidth)
    return axes

def get_mean_value(self, extent=None, step=0):
    if extent is None:
        extent = self.get_extent()
    idx = self.get_extent_idx(extent)
    if type(self.values) == type(list()):
        values = self.values[step][idx]
    else:
        values = self.values[idx]
    return np.mean(values)

def get_values_under_Path(self, path, **kwargs):
    xin = self.x
    yin = self.y
    if len(self.values.shape) == 3:
        zin = self.values[:,0,timestep]
    else:
        zin = self.values
    idx, = np.where(np.logical_and(
        np.logical_and(xin>=np.min(path.vertices[:,0]), xin<=np.max(path.vertices[:,0])),
        np.logical_and(yin>=np.min(path.vertices[:,1]), yin<=np.max(path.vertices[:,1]))))
    return griddata((xin[idx], yin[idx]), zin[idx], (path.vertices[:,0],path.vertices[:,1]))


def interpolate_DEM(self, tile, **kwargs):
    method = kwargs.pop("method", "FVM")
    channel_polygons = kwargs.pop("channel_polygons", None)
    bar_polygons = kwargs.pop("bar_polygons", None)
    tile_extent = tile.get_extent()
    idxs_of_points_in_tile = self.get_extent_idx(tile_extent, tile.epsg)
    xyz = tile.get_xyz(epsg=self.epsg)
    for idx in idxs_of_points_in_tile:
        path = self.get_finite_volume_Path(idx)
        _idx, = np.where(np.logical_and(
                    np.logical_and(xyz[:,0]>=np.min(path.vertices[:,0]), xyz[:,0]<=np.max(path.vertices[:,0])),
                    np.logical_and(xyz[:,1]>=np.min(path.vertices[:,1]), xyz[:,1]<=np.max(path.vertices[:,1]))))
        values = xyz[_idx,:]
        if method=="FVM":
            _values = values[np.where(path.contains_points(values[:,0:2])),2]
            self.values[idx] = np.mean(_values)
        else:
            self.values[idx] = griddata((xyz[_idx,0], xyz[_idx,1]), xyz[_idx,2], (self.x[idx], self.y[idx]), method=method)
        # path.contains_points() takes time, so we try to recycle the _values variable with try/except
        if channel_polygons is not None:
            for channel_polygon in channel_polygons:
                if channel_polygon.contains_point((self.x[idx], self.y[idx])):#path.intersects_path(channel_polygon):
                    try: 
                        _values
                    except:
                        _values = values[np.where(path.contains_points(values[:,0:2])),2]
                    self.values[idx] = np.min(_values)
    
        if bar_polygons is not None:
            for bar_polygon in bar_polygons:
                if bar_polygon.contains_point((self.x[idx], self.y[idx])):#path.intersects_path(bar):
                    try:
                        _values
                    except:
                        _values = values[np.where(path.contains_points(values[:,0:2])),2]
                    self.values[idx] = np.max(_values)

def build_outer_polygon(self):

    boundary_list = list()

    if self.ocean_boundaries is not None:
        for i in range(len(self.ocean_boundaries)):
            vertices = list()
            codes = [Path.MOVETO]
            for j in range(self.ocean_boundaries[i].shape[0]):
                idx = self.ocean_boundaries[i][j]
                vertices.append((self.x[idx], self.y[idx]))
                codes.append(Path.LINETO)
            path = Path(vertices, codes[:-1])
            boundary_list.append(path)

    
    if self.land_boundaries is not None:
       
        for i in range(len(self.land_boundaries)):
            vertices = list()
            codes = [Path.MOVETO]
            for j in range(self.land_boundaries[i][0].shape[0]):
                idx = self.land_boundaries[i][0][j]
                vertices.append((self.x[idx], self.y[idx]))
                codes.append(Path.LINETO)
            path = Path(vertices, codes[:-1])
            boundary_list.append(path)

    if self.inflow_boundaries is not None:
        for i in range(len(self.inflow_boundaries)):
            vertices = list()
            codes = [Path.MOVETO]
            for j in range(self.inflow_boundaries[i][0].shape[0]):
                idx = self.inflow_boundaries[i][0][j]
                vertices.append((self.x[idx], self.y[idx]))
                codes.append(Path.LINETO)
            path = Path(vertices, codes[:-1])
            boundary_list.append(path)
        
    if self.outflow_boundaries is not None:
        for i in range(len(self.outflow_boundaries)):
            vertices = list()
            codes = [Path.MOVETO]
            for j in range(self.outflow_boundaries[i][0].shape[0]):
                idx = self.outflow_boundaries[i][0][j]
                vertices.append((self.x[idx], self.y[idx]))
                codes.append(Path.LINETO)
            path = Path(vertices, codes[:-1])
            boundary_list.append(path)
    
    # if self.weirBoundaries is not None:
        # vertices=list()
        # for boundary in self.weirBoundaries:
            # vertices.append(np.vstack(((self.x[boundary['front_face']], self.y[boundary['front_face']]))).T)
            # vertices.append(np.vstack(((self.x[boundary['back_face']], self.y[boundary['back_face']]))).T)
        # print(vertices)
            
            # vertices=[self.]
            # for face in boundary['front_face']
                
        
            # boundary_list.append(path)
    
    try:
        ordered_list=[boundary_list.pop()]
    except:
        return
        
    b1 = ordered_list[-1]
    b1_bottom_lon =  b1.vertices[-1,0]
    b1_bottom_lat =  b1.vertices[-1,1]

    while boundary_list:
        diff = list()    
        for boundary in boundary_list:
            eudiff = np.sqrt((b1_bottom_lon-boundary.vertices[0,0])**2. + (b1_bottom_lat-boundary.vertices[0,1])**2.)
            diff.append(eudiff)
        # append the boundary with the minimum euclidean difference, which we assume to be next in the sequence
        ordered_list.append(boundary_list.pop(diff.index(min(diff))))
        # get the new boundary limits to find the next shape in the sequence
        b1 = ordered_list[-1]
        b1_bottom_lon = b1.vertices[-1,0]
        b1_bottom_lat = b1.vertices[-1,1]

    final_list = list()
    for path in ordered_list:
        for i in range(len(path.vertices)):
            final_list.append((path.vertices[i,0],path.vertices[i,1]))

    codes = len(final_list) * [Path.LINETO]
    codes[0] = Path.MOVETO
    final_list.append(final_list[0])
    codes.append(Path.CLOSEPOLY)
    return Path(final_list, codes)

def build_inner_polygons(self, epsg=None):
    if self.inner_boundaries is None:
        return
    if epsg is None:
        epsg = self.epsg

    if epsg != self.epsg:
        self_proj = pyproj.Proj(init='epsg:{}'.format(self.epsg))
        target_proj = pyproj.Proj(init='epsg:{}'.format(epsg))
    innerBoundaries = list()
    for i in range(len(self.inner_boundaries)):
        idxs =self.inner_boundaries[i][0]
        vertices = list()
        for idx in idxs:
            vertices.append((self.x[idx], self.y[idx]))
        if epsg != self.epsg:
            x = [x for x,y in vertices]
            y = [y for x,y in vertices]
            x, y = pyproj.transform(self_proj, target_proj, x, y)
            vertices = list(zip(x,y))
        innerBoundaries.append(Path(vertices, closed=True))
        # Uncomment to see a plot of the geometries as they are generated.
        # import matplotlib.pyplot as plt
        # import matplotlib.patches as patches
        # fig = plt.figure()
        # ax = fig.add_subplot(111)
        # patch = patches.PathPatch(path, facecolor='orange', lw=2)
        # ax.add_patch(patch)
        # ax.axis('scaled')
        # plt.show()
    return innerBoundaries

def get_land_boundaries(self, epsg=None):
    if self.land_boundaries is not None:
        if epsg is None:
            epsg = self.epsg
        if epsg != self.epsg:
            self_proj = pyproj.Proj(init='epsg:{}'.format(self.epsg))
            target_proj = pyproj.Proj(init='epsg:{}'.format(epsg))
        boundary_list=list()
        for i in range(len(self.land_boundaries)):
            vertices = list()
            for j in range(self.land_boundaries[i][0].shape[0]):
                idx = self.land_boundaries[i][0][j]
                vertices.append((self.x[idx], self.y[idx]))
            if epsg != self.epsg:
                x = [x for x,y in vertices]
                y = [y for x,y in vertices]
                x, y = pyproj.transform(self_proj, target_proj, x, y)
                vertices = list(zip(x,y))
            boundary_list.append(Path(vertices, closed=True))
        return boundary_list

def plot_outerBoundary(self, extent=None, axes=None, **kwargs):
    axes, idx = _plotters._init_fig(self, axes, extent, title)
    patch = patches.PathPatch(self.outerBoundary, **kwargs)
    axes.add_patch(patch)


def _init_TPXO(self):
    pass