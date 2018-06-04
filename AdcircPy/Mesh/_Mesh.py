import numpy as np
from matplotlib.path import Path
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from haversine import haversine
from AdcircPy.Surface import _fig
from AdcircPy import Mesh
from AdcircPy.Mesh import _fort14

def init_from_fort14(fort14, datum='MSL', epsg=4326):
    return Mesh.Mesh(**_fort14.parse_fort14(fort14))

def plot_bathy(self, extent=None, axes=None, title=None, total_colors=256, cbar_label=r'elevation [m]', **kwargs):
    axes, idx = _fig.init_fig(self, axes, extent, title)
    vmin = kwargs.pop("vmin", np.min(self.values[idx]))
    vmax = kwargs.pop("vmax", np.max(self.values[idx]))
    if vmax < 0.:

        levels = kwargs.pop("levels", np.linspace(vmin, vmax, total_colors))
        mlevel = np.mean(self.values[idx])
        cmap   = kwargs.pop("cmap", plt.cm.seismic)
        col_val = 0.5      
    else:
        mlevel = 0.
        wet_count = int(np.floor(total_colors * \
                (float((self.values[idx] < 0.).sum())/self.values[idx].size)))
        col_val = float(wet_count)/float(total_colors)
        dry_count = total_colors - wet_count
        colors_undersea = plt.cm.bwr(np.linspace(1., 0., wet_count))
        colors_land = plt.cm.terrain(np.linspace(0.25, 1., dry_count))
        colors = np.vstack((colors_undersea, colors_land))
        cmap = LinearSegmentedColormap.from_list('cut_terrain', colors)
        wlevels = np.linspace(vmin, mlevel, wet_count, endpoint=False)
        dlevels = np.linspace(mlevel, vmax, dry_count)
        levels = np.hstack((wlevels, dlevels))
    norm = fig.FixPointNormalize(sealevel=mlevel, vmax=vmax, vmin=vmin, col_val=col_val)
    axes.tricontourf(self.x, self.y, self.elements, self.values,
                    levels=levels, cmap=cmap, norm=norm, extend='both', **kwargs)
    cbar = fig._init_colorbar(axes, cmap, vmin, vmax)
    cbar.set_label(cbar_label)
    cbar.set_ticks([vmin, vmin + col_val *(vmax-vmin), vmax])
    cbar.set_ticklabels([np.around(vmin, 2), mlevel, np.around(vmax, 2)])
    return axes

def plot_shoreline(self, extent=None, axes=None, title=None, color='black', linewidth=0.5):
    axes, idx = _fig.init_fig(self, axes, extent, title)
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
    channel_polygons = kwargs.pop("channel_polygons", None)
    bar_polygons = kwargs.pop("bar_polygons", None)
    tile_extent = tile.get_extent()
    idxs_of_points_in_tile = self.get_extent_idx(extent=tile_extent)
    xyz = tile.get_xyz(epsg=self.epsg)
   
    for i, obs_point_idx in enumerate(idxs_of_points_in_tile):

        adjacent_elements = self.get_elements_surrounding_node(obs_point_idx)
        midpoints = list()
        centroids = list()
        _adjacent_elements=list(adjacent_elements)

        # Iterate over surrounding elements
        for j, element in enumerate(_adjacent_elements):
            _element = list(element)
            _element.append(_element[0])
            _midpoints = list()
            
            # Calculate element midpoints and centroid
            for i, point_idx in enumerate(element):         
                dx =  self.x[_element[i+1]] - self.x[_element[i]]
                dy =  self.y[_element[i+1]] - self.y[_element[i]]
                if np.isin(obs_point_idx, [_element[i], _element[i+1]]):
                    midpoints.append((self.x[_element[i]] + 0.5*dx, self.y[_element[i]] + 0.5*dy))
                _midpoints.append((self.x[_element[i]] + 0.5*dx, self.y[_element[i]] + 0.5*dy))
            _midpoints = np.array(_midpoints)
            centroids.append((np.mean(_midpoints[:,0]),np.mean(_midpoints[:,1])))
            _element = list(element)
            while _element[0] != obs_point_idx:
                _element = list(np.roll(_element,1))
            _adjacent_elements[j] = _element
        _adjacent_elements = np.array(_adjacent_elements)
        lateral_indexes, counts = np.unique(_adjacent_elements, return_counts=True)
        vertices = list()
        for centroid in centroids:
            vertices.append(centroid)
        for midpoint in midpoints:
            vertices.append(midpoint)    
        ordered_vertices = list()
        if 1 in counts: # means that this node is at a boundary and need to be included in the polygon.
            ordered_vertices.append((self.x[obs_point_idx], self.y[obs_point_idx]))
            index = list(np.delete(lateral_indexes, np.where(counts>1))).pop()
            x = self.x[obs_point_idx]
            y = self.y[obs_point_idx]
            dx = self.x[index] - self.x[obs_point_idx]
            dy = self.y[index] - self.y[obs_point_idx]
            ordered_vertices.append(vertices.pop(vertices.index((x+0.5*dx,y+0.5*dy))))
        else: # means that this is node is fully surrounded by elements so the obs point is not included in polygon.
            ordered_vertices.append(vertices.pop())

        lon = ordered_vertices[-1][0]
        lat = ordered_vertices[-1][1]
        while vertices:
            diff = list()
            for vertex in vertices:
                _diff = haversine((vertex[1], vertex[0]), (lat, lon))
                diff.append(_diff)
            ordered_vertices.append(vertices.pop(diff.index(min(diff))))
            lon = ordered_vertices[-1][0]
            lat = ordered_vertices[-1][1]
        ordered_vertices.append(ordered_vertices[0])
        path = Path(ordered_vertices, closed=True)
        points = np.asarray(ordered_vertices)
        idx, = np.where(np.logical_and(
                            np.logical_and(xyz[:,0]>=np.min(points[:,0]), xyz[:,0]<=np.max(points[:,0])),
                            np.logical_and(xyz[:,1]>=np.min(points[:,1]), xyz[:,1]<=np.max(points[:,1]))))
        values = xyz[idx,:]
        values = values[np.where(path.contains_points(values[:,0:2])),2]
        self.values[obs_point_idx] = np.mean(values)
        
        if channel_polygons is not None:
            for channel_polygon in channel_polygons:
                if channel_polygon.contains_point((self.x[obs_point_idx], self.y[obs_point_idx])):#path.intersects_path(channel_polygon):
                    self.values[obs_point_idx] = np.min(values)
        
        if bar_polygons is not None:
            for bar_polygon in bar_polygons:
                if bar_polygon.contains_point((self.x[obs_point_idx], self.y[obs_point_idx])):#path.intersects_path(bar):
                    self.values[obs_point_idx] = np.max(values)     
        
        
        # # Uncomment here to make a plot of the calculations we just made.
        # # Put under statement "if 1 in counts:" if you only want to see boundary nodes.
        # # if 1 in counts:
        # import matplotlib.patches as patches
        # fig = plt.figure()
        # ax  = fig.add_subplot(111)
        # in_polygon  = path.contains_points(xyz[:,0:2])
        # in_polygon  = in_polygon.reshape(tile.values.shape)
        # x,y = np.meshgrid(tile.x, tile.y)
        # z = np.flipud(tile.values)
        # z = np.ma.masked_where(~in_polygon, z)
        # ax.pcolormesh(x, y, np.flipud(z),cmap='jet')
        # lon = list()
        # lat = list()
        # for element in _adjacent_elements:
            # vertices = list()
            # for index in element:
                # lon.append(self.x[index])
                # lat.append(self.y[index])
                # vertices.append((self.x[index],self.y[index]))
            # vertices.append(vertices[0])
            # _element_as_path = Path(vertices, closed=True)
            # ax.add_patch(patches.PathPatch(_element_as_path, facecolor='none', lw=2, edgecolor='black'))
        # ax.add_patch(patches.PathPatch(path, facecolor='none', lw=1, edgecolor='red'))
        # ax.axis('scaled')
        # ax.axis([np.min(lon), np.max(lon), np.min(lat), np.max(lat)])
        # plt.show()

        
        
        
    
    # # Correction for Weir and Culvert boundaries. 
    # temporarily converts non-MSL grids to MSL.
    # if self.datum != 'MSL':
    #     self.values = self.values + self.datum_offset
        
    # for boundary in self.weirBoundaries:
    #     for i, idx in enumerate(boundary[0]):
    #         if idx in idxs_of_points_in_tile and \
    #             self.values[idx] > boundary[2][i]:
    #             self.values[idx] = boundary[2][i]
    #     for i, idx in enumerate(boundary[1]):
    #         if idx in idxs_of_points_in_tile and \
    #             self.values[idx] > boundary[2][i]:
    #             self.values[idx] = boundary[2][i]

    # # Need an example to test culvert implementation.
    # for boundary in self.culvertBoundaries:
    #     for i, idx in enumerate(boundary[0]):
    #         if self.values[idx] > boundary[2][i]:
    #             self.values[idx] = boundary[2][i]
    #     for i, idx in enumerate(boundary[1]):
    #         if self.values[idx] > boundary[2][i]:
    #             self.values[idx] = boundary[2][i]
    
    # changes back from MSL to original datum
    # if self.datum != 'MSL':
    #     self.values = self.values - self.datum_offset


    # def _interpolate_DEM(self, tile, radius=0.001):
    # """
    # This function is a prototype and should not be used. It is only in the
    # source code because it contains some techniques that could potentially
    # be useful for other uses in the future and is provided merely as reference.
    # This function will be entirely removed in future versions.
    # """
    # xyz = tile.get_xyz()
    
    # xs  = xyz[:,0]
    # ys  = xyz[:,1]
    # zs  = xyz[:,2]
    # tile_extent = tile.get_extent()
    # tri_paths   = self.get_elements_as_paths(extent=tile_extent)
    # idx         = self.get_extent_idx(extent=tile_extent)
    # self.values[idx] = -99999.0
    
    
    # # Begin channel detection algorithm
    # structure = ndimage.generate_binary_structure(2,2)
    # for i, triangle in enumerate(tri_paths):
        # in_tri = triangle.contains_points(xyz[:,:2])
        # total_idx, = np.where(in_tri)
        # wet_idx, = np.where(np.logical_and(in_tri, zs<0.0))
        # dry_idx, = np.where(np.logical_and(in_tri, zs>=0.0))

        # # if there are both wet and dry points inside the triangle
        # if len(wet_idx)>0 and len(dry_idx)>0 and \
            # float(zs[wet_idx].size)/float(zs[total_idx].size)>=0.5:
            
            # dry_bool_lg = np.logical_and(
                            # triangle.contains_points(xyz[:,:2], radius=radius),
                            # zs >= 0.0)
            
            # dry_idx_lg = np.where(dry_bool_lg.reshape(tile.values.shape))
            # image = np.zeros(tile.values.shape)
            # image[dry_idx_lg] = 1.0
            # blob_img, num_blobs = ndimage.label(image, structure=structure)
            # node_neighbor_idx, = np.where(dry_bool_lg)
            # node_ids = griddata((xs[node_neighbor_idx], ys[node_neighbor_idx]),
                                # blob_img.reshape(xs.shape)[node_neighbor_idx],
                            # (triangle.vertices[:-1,0], triangle.vertices[:-1,1]),
                                # method='nearest')
            # dry_nodes = filter(lambda x: x!=0., node_ids)
            # unique_ids = len(np.unique(dry_nodes))
            
            # if unique_ids>=2 and np.mean(zs[wet_idx])<-1.:
                
                # vertices_idx, = np.where(np.logical_and(
                    # np.isin(self.x, triangle.vertices[:-1,0]),
                    # np.isin(self.y, triangle.vertices[:-1,1])))
                # self.values[vertices_idx] = np.mean(zs[wet_idx])


    # idx, = np.where(np.isin(self.values, -99999.0))
    # values = griddata((xs, ys), zs, (self.x[idx], self.y[idx]), method='linear')    

    # for p, index in enumerate(idx):
        # self.values[index] = values[p]
    
    # # temporarily converts non-MSL grids to MSL.
    # if self.datum != 'MSL': self.values = self.values + self.datum_offset
        
        
    # paths = self.get_contours(levels=[0], extent=tile_extent)
    # for path in paths:
        # if 79 in path.codes:
            # idx, = np.where(np.logical_and(
                        # path.contains_points(np.vstack((self.x, self.y)).T),
                        # self.values<0.))
            # self.values[idx]=0.
    
    # # changes back from MSL to original datum
    # if self.datum != 'MSL': self.values = self.values - self.datum_offset
    

