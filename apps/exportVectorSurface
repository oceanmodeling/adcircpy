from __future__ import absolute_import, division, print_function
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def vectorize(self, **kwargs):
    
    levels = kwargs.pop("levels", self.levels)
    if np.ma.is_masked(self.values):
        # trimask = np.any(self.values.mask[self.elements], axis=1)
        point_mask_indices = np.where(self.values.mask)
        trimask = np.any(np.in1d(self.elements,point_mask_indices).reshape(-1,3),axis=1)
        Tri = matplotlib.tri.Triangulation(self.x, self.y, self.elements, trimask)
    else:
        Tri = matplotlib.tri.Triangulation(self.x, self.y, self.elements)
        
    contour = plt.tricontourf(Tri, self.values, levels=levels)
    # plt.show()

    for Level in contour.collections:
        for Path in Level.get_paths():
            # Path.simplify_threshold = 0.0
            innerPolygons = list()
            outerPolygons = list()
            polygons = [x for x in Path.to_polygons() if x.shape[0] >=3]
            print(len(polygons))
            BREAK
            for polygon in Path.to_polygons():
                if np.cross(polygon, np.roll(polygon, -1, axis=0)).sum() >= 0:
                    outerPolygons.append(polygon)
                else:
                    innerPolygons.append(polygon)
                print(len(outerPolygons))
                print(len(innerPolygons))
                BREAK
        # for outerPolygon in outerPolygons:
        #     for innerPolygon in innerPolygons:
        #         path = matplotlib.path.Path(outerPolygon)
        #         PointsInside = path.contains_points(innerPolygon)
                
        #         if  PointsInside.any() == True:
        #             if PointsInside.all() == True:
        #                 print("All points inside")
        #             else:
        #                 print("Only some of them")

            # inner_points = [x[0] for x in innerPolygons]
            # print(inner_points[0])
            
            # inout = pointsInsidePoly(inner_points,out)
            # overall_inout = np.logical_or(overall_inout, inout)
            # out_inner = [g for f, g in enumerate(inner) if inout[f]]
            # poly = Polygon(out, out_inner)
            # pass
        # else:
            # poly = Polygon(out)#! /usr/bin/env python


            def get_contours(self, levels, extent=None, timestep=0):
    if len(self.values.shape) == 3:
        zin = self.values[:,0,timestep]
    else:
        zin = self.values
    fig, ax = plt.subplots()
    cs = ax.tricontour(self.x, self.y, self.elements, zin, levels=levels)
    plt.draw()
    plt.close(fig)
    paths = cs.collections[0].get_paths()
    for i, path in enumerate(list(paths)):
        codes = len(path.vertices) * [Path.LINETO]
        codes[0] = Path.MOVETO
        if len(path.vertices) >= 3 and \
                (path.vertices[0,0] == path.vertices[-1,0] and \
                    path.vertices[0,-1] == path.vertices[-1,-1]):
            codes[-1] = Path.CLOSEPOLY
        paths[i] = Path(path.vertices, codes)
    if extent is not None:
        extentWindow =  Path([(extent[0], extent[2]), (extent[1], extent[2]),
        (extent[1], extent[3]), (extent[0], extent[3]),(extent[0], extent[2])], 
        [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY])
        for i, path in enumerate(list(paths)):
            if not extentWindow.intersects_path(path, filled=False):
                paths.remove(path)
    return np.asarray(paths).flatten()