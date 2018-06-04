from __future__ import absolute_import, division, print_function
import numpy as np
from matplotlib.path import Path
import pyproj

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