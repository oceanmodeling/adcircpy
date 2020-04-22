from matplotlib.tri import Triangulation
import numpy as np
from pyproj import CRS
import matplotlib.pyplot as plt
from collections import defaultdict
from itertools import permutations
from matplotlib.transforms import Bbox
from pyproj import Transformer


class UnstructuredMesh:

    """Graph-like representation of a geographical mesh.

        Arguments x and y are iterables of the same length.

        >>> from adcircpy.mesh import UnstructuredMesh
        >>> x = [0, 1, 0.5]
        >>> y = [0, 0, 0.5]
        >>> mesh = UnstructuredMesh(x, y)
        >>> isinstance(mesh, UnstructuredMesh)
        True

        The elements argument is optional:

        >>> elements = [[0, 1, 2]]
        >>> mesh = UnstructuredMesh(x, y, elements=elements)
        >>> isinstance(mesh, UnstructuredMesh)
        True

    """

    def __init__(
        self,
        vertices,
        elements=None,
        crs=None,
    ):
        """
        :param iterable x: an iterable of x coordinates
        :param iterable y: an iterable of y coordinates
        :param iterable triangles: an iterable of triangle indexing
        :param str crs: a proj4 init string
        """
        self._vertices = vertices
        self._elements = elements
        self._crs = crs

    @classmethod
    def from_triangulation(cls, triangulation, crs=None):
        """
        Alternative constructor.

        :param triangulation: triangulation object
        :type Triangulation: matplotlib.tri.Triangulation instance

        >>> from matplotlib.tri import Triangulation
        >>> tri = Triangulation(mesh.x, mesh.y, mesh.triangles)
        >>> mesh = mesh.from_triangulation(tri)
        """
        return cls(
            np.vstack([triangulation.x, triangulation.y]),
            triangulation.triangles,
            crs)

    def add_attribute(self, name, **properties):
        if self.has_attribute(name):
            raise AttributeError(
                'Non-unique attribute name: '
                + 'Attribute attribute name already exists.')
        else:
            self._attributes[name] = dict()
            self._attributes[name]['values'] = None
            self._attributes[name].update(properties)

    def has_attribute(self, name):
        if name in self._attributes.keys():
            return True
        else:
            return False

    def get_attribute(self, name):
        if not self.has_attribute(name):
            raise AttributeError(f'Attribute {name} not set.')
        return self._attributes[name]

    def get_attribute_values(self, name):
        if not self.has_attribute(name):
            raise AttributeError(f'Attribute {name} not set.')
        return self._attributes[name]['values']

    def get_attribute_properties(self, name):
        if not self.has_attribute(name):
            raise AttributeError(f'Attribute {name} not set.')
        return self._attributes[name]['properties']

    def get_attribute_names(self):
        return list(self._attributes.keys())

    def set_attribute(self, name, values, elements=False, **properties):
        if name not in self.get_attribute_names():
            raise AttributeError(
                f'Cannot set attribute: {name} is not an attribute.')
        values = np.asarray(values)
        assert isinstance(elements, bool)
        if elements:
            assert values.shape[0] == self.elements.shape[0]
        else:
            assert values.shape[0] == self.vertices.shape[0]
        self._attributes[name]['values'] = values
        self._attributes[name].update(properties)

    def remove_attribute(self, name):
        if name in self.get_attribute_names():
            self._attributes.pop(name)
        else:
            raise AttributeError(
                'Cannot remove attribute: attribute does not exist.')

    def transform_to(self, dst_crs):
        dst_crs = CRS.from_user_input(dst_crs)
        if self.srs != dst_crs.srs:
            transformer = Transformer.from_crs(
                self.crs, dst_crs, always_xy=True)
            x, y = transformer.transform(self.x, self.y)
            self._vertices = np.vstack([x, y]).T
            self._crs = dst_crs

    def triplot(
        self,
        ax=None,
        show=False,
        title=None,
        color='k',
        linewidth=0.5,
        **kwargs
    ):
        kwargs.update({
            'color': color,
            'linewidth': linewidth
            })
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        ax.triplot(self.triangulation, **kwargs)
        if title:
            plt.title(title)
        if show:
            ax.axis('scaled')
            plt.show()

    def assign_crs(self, crs):
        self._crs = crs

    def get_x(self, crs=None):
        return self.get_xy(crs)[:, 0]

    def get_y(self, crs=None):
        return self.get_xy(crs)[:, 1]

    def get_xy(self, crs=None):
        if crs is not None:
            crs = CRS.from_user_input(crs)
            transformer = Transformer.from_crs(self.crs, crs, always_xy=True)
            x, y = transformer.transform(self.x, self.y)
            return np.vstack([x, y]).T
        else:
            return np.vstack([self.x, self.y]).T

    def get_xyz(self, crs=None):
        return np.hstack([self.get_xy(crs), self.values])

    def get_extent(self, crs=None):
        xy = self.get_xy(crs)
        return (np.min(xy[:, 0]), np.max(xy[:, 0]),
                np.min(xy[:, 1]), np.max(xy[:, 1]))

    def get_bbox(self, crs=None):
        xmin, xmax, ymin, ymax = self.get_extent(crs)
        return Bbox([[xmin, ymin], [xmax, ymax]])

    @property
    def vertices(self):
        return self._vertices

    @property
    def x(self):
        return self.vertices[:, 0]

    @property
    def y(self):
        return self.vertices[:, 1]

    @property
    def triangulation(self):
        try:
            return self.__triangulation
        except AttributeError:
            self.__triangulation = Triangulation(self.x, self.y, self.elements)
            return self.__triangulation

    @property
    def elements(self):
        return self._elements

    @property
    def node_neighbors(self):
        return self._node_neighbors

    @property
    def proj(self):
        return CRS.from_user_input(self.crs)

    @property
    def crs(self):
        return self._crs

    @property
    def bbox(self):
        return self.get_bbox()

    @property
    def _attributes(self):
        try:
            return self.__attributes
        except AttributeError:
            self.__attributes = dict()
            return self.__attributes

    @property
    def _vertices(self):
        return self.__vertices

    @property
    def _elements(self):
        return self.__elements

    @property
    def _crs(self):
        return self.__crs

    @property
    def _node_neighbors(self):
        try:
            return self.__node_neighbors
        except AttributeError:
            self.__node_neighbors = defaultdict(set)
            for simplex in self.triangulation.triangles:
                for i, j in permutations(simplex, 2):
                    self.__node_neighbors[i].add(j)
            return self.__node_neighbors

    @_vertices.setter
    def _vertices(self, vertices):
        vertices = np.asarray(vertices)
        msg = "ERROR: vertices array must be of shape (N, 2)"
        assert vertices.shape[1] == 2, msg
        self.__vertices = vertices

    @_elements.setter
    def _elements(self, elements):
        if elements is None:
            elements = Triangulation(self.x, self.y).triangles
        self.__elements = np.asarray(elements)

    @_crs.setter
    def _crs(self, crs):
        if not isinstance(crs, (type(None), CRS)):
            crs = CRS.from_user_input(crs)
        self.__crs = crs


if __name__ == "__main__":
    import doctest
    x = [0, 1, 0.5]
    y = [0, 0, 0.5]
    mesh = UnstructuredMesh(x, y)
    doctest.testmod(extraglobs={'mesh': mesh})


























































# class UnstructuredMesh:

#     def __init__(self, vertices, elements, values=None, SpatialReference=None):
#         super(UnstructuredMesh, self).__init__()
#         self._vertices = vertices
#         self._elements = elements
#         self._values = values
#         self._SpatialReference = SpatialReference
#         self.__attributes = dict()



#     def interpolate(self, Dataset):
#         assert isinstance(Dataset, gdal.Dataset)
#         if not self.SpatialReference.IsSame(
#                     gdal_tools.get_SpatialReference(Dataset)):
#             Dataset = gdal_tools.Warp(Dataset, dstSRS=self.SpatialReference)
#         x, y, z = gdal_tools.get_arrays(Dataset)
#         bbox = gdal_tools.get_Bbox(Dataset)
#         f = RectBivariateSpline(x, y, z.T, bbox=[bbox.xmin, bbox.xmax,
#                                                  bbox.ymin, bbox.ymax])
#         idxs = np.where(np.logical_and(
#                             np.logical_and(
#                                 bbox.xmin <= self.vertices[:, 0],
#                                 bbox.xmax >= self.vertices[:, 0]),
#                             np.logical_and(
#                                 bbox.ymin <= self.vertices[:, 1],
#                                 bbox.ymax >= self.vertices[:, 1])))[0]
#         values = f.ev(self.vertices[idxs, 0], self.vertices[idxs, 1])
#         new_values = self.values.copy()
#         for i, idx in enumerate(idxs):
#             new_values[idx] = values[i]
#         self.values = new_values

#     def has_invalid(self):
#         return np.any(np.isnan(self.values))

#     def fix_invalid(self, method='nearest'):
#         if self.has_invalid():
#             if method == 'nearest':
#                 idx = np.where(~np.isnan(self.values))
#                 _idx = np.where(np.isnan(self.values))
#                 values = griddata(
#                     (self.x[idx], self.y[idx]), self.values[idx],
#                     (self.x[_idx], self.y[_idx]), method='nearest')
#                 new_values = self.values.copy()
#                 for i, idx in enumerate(_idx):
#                     new_values[idx] = values[i]
#                 self._values = new_values
#                 return self.values
#             else:
#                 raise NotImplementedError

#     def make_plot(self, show=False, levels=256):
#         z = np.ma.masked_invalid(self.values)
#         vmin, vmax = z.min(), z.max()
#         z = z.filled(fill_value=-99999.)
#         if isinstance(levels, int):
#             levels = np.linspace(vmin, vmax, levels)
#         plt.tricontourf(self.mpl_tri, z, levels=levels)
#         plt.gca().axis('scaled')
#         if show:
#             plt.show()
#         plt.gca().axis('scaled')
#         return plt.gca()



#     @property
#     def vertices(self):
#         return self.__vertices

#     @property
#     def elements(self):
#         return self.__elements

#     @property
#     def values(self):
#         return self.__values

#     @property
#     def x(self):
#         return self.vertices[:, 0]

#     @property
#     def y(self):
#         return self.vertices[:, 1]

#     @property
#     def xy(self):
#         return self.vertices

#     @property
#     def xyz(self):
#         return self.get_xyz()

#     @property
#     def SpatialReference(self):
#         return self._get_spatial_reference()

#     @property
#     def planar_straight_line_graph(self):
#         """
#         Slow alogirthm, open for suggestions...
#         """
#         try:
#             return self.__planar_straight_line_graph
#         except AttributeError:
#             _pslg = _PlanarStraightLineGraph(self.SpatialReference)
#             idxs = np.vstack(list(np.where(self.mpl_tri.neighbors == -1))).T
#             unique_edges = list()
#             for i, j in idxs:
#                 unique_edges.append((self.mpl_tri.triangles[i, j],
#                                      self.mpl_tri.triangles[i, (j+1) % 3]))
#             unique_edges = np.asarray(unique_edges)
#             ring_collection = list()
#             initial_idx = 0
#             for i in range(1, len(unique_edges)-1):
#                 if unique_edges[i-1, 1] != unique_edges[i, 0]:
#                     try:
#                         idx = np.where(
#                             unique_edges[i-1, 1] == unique_edges[i:, 0])[0][0]
#                         unique_edges[[i, idx+i]] = unique_edges[[idx+i, i]]
#                     except IndexError:
#                         ring_collection.append(unique_edges[initial_idx:i, 0])
#                         initial_idx = i
#                         continue
#             if len(ring_collection) == 0:
#                 ring_collection.append(unique_edges[initial_idx:i, 0])
#             print(
#                 [Path(np.hstack([self.x[ring], self.y[ring]]).T, closed=True)
#                  for ring in ring_collection])
#             BREAKME


#             # linear_ring_collection = list()
#             # for ring in ring_collection:
#             #     # plt.plot()
#             #     _geom = ogr.Geometry(ogr.wkbLinearRing)
#             #     _geom.AssignSpatialReference(self.SpatialReference)
#             #     for idx in ring:
#             #         _geom.AddPoint(self.x[idx], self.y[idx])
#             #     _geom.CloseRings()
#             #     _polygon = ogr.Geometry(ogr.wkbPolygon)
#             #     _polygon.AddGeometry(_geom)
#             #     linear_ring_collection.append(_polygon)
#             # idx = np.where(
#             #     np.max([geom.GetArea() for geom in linear_ring_collection]))[0]
#             # linear_ring = linear_ring_collection.pop(idx[0])
#             # # contains = [linear_ring.Contains(geom)
#             # #             for geom in linear_ring_collection]
#             # contains = [geom.Within(linear_ring)
#             #             for geom in linear_ring_collection]
#             # print(contains)
#             # BREAKME

#                 # print(_polygon.IsValid())
#                 # _geom.CloseRings()
#                 # print(_geom.IsValid())
#                 # BREAKME
#                 # vertices = np.asarray([(x, y) for x, y, _ in _geom.GetPoints()])
#                 # plt.plot(vertices[:, 0], vertices[:, 1])
#                 # plt.show()

#                 # _geom.AssignSpatialReference(self.SpatialReference)
#                 # print(_geom.IsValid())
#                 # BREAKME
#                 # linear_ring_collection.append(_geom)
#             # print([geom.Polygonize() for geom in linear_ring_collection])
#             # print([geom.IsRing() for geom in linear_ring_collection])
#             # BREAKME
#             # idx = np.where(
#             #     np.max([geom.GetArea() for geom in linear_ring_collection]))[0]
#             # linear_ring = linear_ring_collection.pop(idx[0])
#             # contains = [linear_ring.Contains(geom)
#             #             for geom in linear_ring_collection]
#             # print(contains)
#             # BREAKME

#         # print(linear_ring)

#         # while len(linear_ring_collection) > 0:
#         #     linear_ring = linear_ring_collection.pop()
#         #     for remaining_rings in linear_ring_collection:


                



#         # print(geom_collection)
#         # BREAKME

#         # ring_collection = []
#         # while len(ring_collection) > 0:
#         #     ring = ring_collection.pop()
#         #     for test_ring in ring_collection:
#         #         if ring.contains(test_ring)
#         # multipolygon = ogr.Geometry(ogr.wkbMultiPolygon)
#         # for polygon in polygons:
#         #     multipolygon.AddGeometry(polygon)
#         # #  -------------------
#         # geom_collection = list()
#         # for ring in ring_collection:
#         #     _geom = ogr.Geometry(ogr.wkbLinearRing)
#         #     _geom.AssignSpatialReference(self.SpatialReference)
#         #     for idx in ring:
#         #         _geom.AddPoint_2D(self.x[idx[0]], self.y[idx[0]])
#         #     geom_collection.append(_geom)
#         # lengths = [_geom.Length() for _geom in geom_collection]
#         # outer_edges = ring_collection.pop(
#         #         np.where(np.max(lengths) == lengths)[0][0])
#         # inner_edges = ring_collection
#         # outer_vertices = self.vertices[outer_edges[:, 0]]
#         # outer_vertices = np.vstack([outer_vertices, outer_vertices[0, :]])
#         # inner_vertices = [self.vertices[ring[:, 0]] for ring in inner_edges]
#         # inner_vertices = [np.vstack([vertices, vertices[0, :]])
#         #                   for vertices in inner_vertices]
#         # return _PlanarStraightLineGraph(
#         #         self.SpatialReference, outer_vertices, *inner_vertices,
#         #         outer_edges=outer_edges, inner_edges=inner_edges)
#             return self.__planar_straight_line_graph

#     # def make_plot(self):
#     #     plt.plot(self.outer_vertices[:, 0], self.outer_vertices[:, 1])
#     #     for vertices in self.inner_vertices:
#     #         plt.plot(vertices[:, 0], vertices[:, 1])
#     #     plt.show()

#     @property
#     def mpl_tri(self):
#         # Since vertices, elements are constants, Triangulation cannot change.
#         try:
#             return self.__mpl_tri
#         except AttributeError:
#             if not np.ma.is_masked(self.values):
#                 self.__mpl_tri = Triangulation(self.x, self.y, self.elements)
#                 return self.__mpl_tri
#             else:
#                 self.__mpl_tri = Triangulation(
#                     self.x, self.y, self.elements,
#                     mask=np.any(self.values.mask[self.elements], axis=1))
#                 return self.__mpl_tri

#     @property
#     def ndim(self):
#         return 2

#     @property
#     def num_elements(self):
#         return self.elements.shape[0]

#     @property
#     def num_nodes(self):
#         return self.vertices.shape[0]

#     @SpatialReference.setter
#     def SpatialReference(self, SpatialReference):
#         assert isinstance(SpatialReference, (int, osr.SpatialReference)), \
#             "Input must be a EPSG code or osr.SpatialReference instance."
#         msg = "Mesh must have a spatial reference assigned before "
#         msg += "transformation can occur."
#         assert self._get_spatial_reference() is not None, msg
#         vertices = self.transform_vertices(
#                 self.__vertices, self._get_spatial_reference(),
#                 SpatialReference)
#         if SpatialReference.IsGeographic() and self._get_spatial_reference().IsProjected():
#             vertices = np.fliplr(vertices)
#         self._vertices = vertices
#         self._SpatialReference = SpatialReference

#     @property
#     def _vertices(self):
#         return self.__vertices

#     @property
#     def _elements(self):
#         return self.__elements

#     @property
#     def _values(self):
#         return self.__values

#     @property
#     def _SpatialReference(self):
#         return self._get_spatial_reference()

#     @_vertices.setter
#     def _vertices(self, vertices):
#         vertices = np.asarray(vertices)
#         assert vertices.shape[1] == self.ndim
#         self.__vertices = vertices

#     @_elements.setter
#     def _elements(self, elements):
#         elements = np.asarray(elements)
#         assert elements.shape[1] == 3
#         self.__elements = elements

#     @_values.setter
#     def _values(self, values):
#         if values is None:
#             values = np.full((self.vertices.shape[0],), np.nan)
#         if not np.ma.is_masked(values):
#             values = np.asarray(values)
#         assert values.shape[0] == self.vertices.shape[0]
#         self.__values = values

#     @_SpatialReference.setter
#     def _SpatialReference(self, SpatialReference):
#         if SpatialReference is not None:
#             self._set_spatial_reference(SpatialReference)
#         else:
#             self._clear_spatial_reference()

        #     # Assume EOF if NOPE is empty.
        #     try:
        #         NOPE = int(f.readline().split()[0])
        #     except IndexError:
        #         return f
        #     # For now, let NOPE=-1 mean a self closing mesh
        #     # reassigning NOPE to 0 until further implementation is applied.
        #     if NOPE == -1:
        #         NOPE = 0
        #     _NOPE = len([])
        #     f.readline()  # Number of total open ocean nodes. Not used.
        #     while _NOPE < NOPE:
        #         grd['ocean_boundaries'].append({'indexes': list()})
        #         NETA = int(f.readline().split()[0])
        #         _NETA = len([])
        #         while _NETA < NETA:
        #             grd['ocean_boundaries'][_NOPE]['indexes'].append(
        #                                         int(f.readline().split()[0])-1)
        #             _NETA += 1
        #         _NOPE += 1
        #     NBOU = int(f.readline().split()[0])
        #     _NBOU = len([])
        #     f.readline()
        #     while _NBOU < NBOU:
        #         NVELL, IBTYPE = map(int, f.readline().split()[:2])
        #         _NVELL = 0
        #         if IBTYPE in [0, 10, 20]:
        #             grd['land_boundaries'].append(
        #                 {'ibtype': IBTYPE,
        #                  'indexes': list()})
        #         elif IBTYPE in [1, 11, 21]:
        #             grd['inner_boundaries'].append(
        #                 {'ibtype': IBTYPE,
        #                  'indexes': list()})
        #         elif IBTYPE in [2, 12, 22, 102, 122]:
        #             grd['inflow_boundaries'].append(
        #                 {'ibtype': IBTYPE,
        #                  'indexes': list()})
        #         elif IBTYPE in [3, 13, 23]:
        #             grd['outflow_boundaries'].append(
        #                 {'ibtype': IBTYPE,
        #                  'indexes': list(),
        #                  'barrier_heights': list(),
        #                  'supercritical_flow_coefficients': list()})

        #         elif IBTYPE in [4, 24]:
        #             grd['weir_boundaries'].append(
        #                 {'ibtype': IBTYPE,
        #                  'front_face_indexes': list(),
        #                  'back_face_indexes': list(),
        #                  'barrier_heights': list(),
        #                  'subcritical_flow_coefficients': list(),
        #                  'supercritical_flow_coefficients': list()})
        #         elif IBTYPE in [5, 25]:
        #             grd['culvert_boundaries'].append(
        #                 {'ibtype': IBTYPE,
        #                  'front_face_indexes': list(),
        #                  'back_face_indexes': list(),
        #                  'barrier_heights': list(),
        #                  'subcritical_flow_coefficients': list(),
        #                  'supercritical_flow_coefficients': list(),
        #                  'cross_barrier_pipe_heights': list(),
        #                  'friction_factors': list(),
        #                  'pipe_diameters': list()})
        #         else:
        #             raise Exception('IBTYPE={} '.format(IBTYPE)
        #                             + 'found in fort.14 not recongnized. ')
        #         while _NVELL < NVELL:
        #             line = f.readline().split()
        #             if IBTYPE in [0, 10, 20]:
        #                 grd['land_boundaries'][-1][
        #                     'indexes'].append(int(line[0])-1)
        #             elif IBTYPE in [1, 11, 21]:
        #                 grd['inner_boundaries'][-1][
        #                     'indexes'].append(int(line[0])-1)
        #             elif IBTYPE in [3, 13, 23]:
        #                 grd['outflow_boundaries'][-1][
        #                     'indexes'].append(int(line[0])-1)
        #                 grd['outflow_boundaries'][-1][
        #                     'external_barrier_heights'].append(float(line[1]))
        #                 grd['outflow_boundaries'][-1][
        #                     'supercritical_flow_coefficients'].append(
        #                         float(line[2]))
        #             elif IBTYPE in [2, 12, 22, 102, 122]:
        #                 grd['iflowBoundaries'][-1][
        #                     'indexes'].append(int(line[0])-1)
        #             elif IBTYPE in [4, 24]:
        #                 grd['weir_boundaries'][-1][
        #                     'front_face_indexes'].append(int(line[0])-1)
        #                 grd['weir_boundaries'][-1][
        #                     'back_face_indexes'].append(int(line[1])-1)
        #                 grd['weir_boundaries'][-1][
        #                     'barrier_heights'].append(float(line[2]))
        #                 grd['weir_boundaries'][-1][
        #                     'subcritical_flow_coefficients'].append(
        #                         float(line[3]))
        #                 grd['weir_boundaries'][-1][
        #                     'supercritical_flow_coefficients'].append(
        #                         float(line[4]))
        #             elif IBTYPE in [5, 25]:
        #                 grd['culvert_boundaries'][-1][
        #                     'front_face_indexes'].append(int(line[0])-1)
        #                 grd['culvert_boundaries'][-1][
        #                     'back_face_indexes'].append(int(line[1])-1)
        #                 grd['culvert_boundaries'][-1][
        #                     'barrier_heights'].append(float(line[2]))
        #                 grd['culvert_boundaries'][-1][
        #                     'subcritical_flow_coefficients'].append(
        #                         float(line[3]))
        #                 grd['culvert_boundaries'][-1][
        #                     'supercritical_flow_coefficients'].append(
        #                         float(line[4]))
        #                 grd['culvert_boundaries'][-1][
        #                     'friction_factors'].append(float(line[5]))
        #                 grd['culvert_boundaries'][-1][
        #                     'pipe_diameters'].append(float(line[6]))
        #             else:
        #                 Exception("Duck-typing error. "
        #                           + "This exception should be unreachable.")
        #             _NVELL += 1
        #         _NBOU += 1
        # return grd