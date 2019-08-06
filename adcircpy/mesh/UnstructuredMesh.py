import numpy as np
from osgeo import osr, gdal, ogr
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.path import Path
from scipy.interpolate import RectBivariateSpline, griddata
from adcircpy.lib import gdal_tools
from adcircpy.lib._SpatialReference import _SpatialReference
from adcircpy.mesh import PlanarStraightLineGraph as _PlanarStraightLineGraph


class UnstructuredMesh(_SpatialReference):

    def __init__(self, vertices, elements, values=None, SpatialReference=None):
        super(UnstructuredMesh, self).__init__()
        self._vertices = vertices
        self._elements = elements
        self._values = values
        self._SpatialReference = SpatialReference
        self.__attributes = dict()

    def get_x(self, SpatialReference=None):
        """ """
        return self.get_xy(SpatialReference)[:, 0]

    def get_y(self, SpatialReference=None):
        """ """
        return self.get_xy(SpatialReference)[:, 1]

    def get_xy(self, SpatialReference=None):
        return self.transform_vertices(self.xy, self.SpatialReference,
                                       SpatialReference)

    def get_xyz(self, SpatialReference=None):
        return np.hstack([self.get_xy(SpatialReference), self.values])

    def get_extent(self, SpatialReference=None):
        xy = self.get_xy(SpatialReference)
        return (np.min(xy[:, 0]), np.max(xy[:, 0]),
                np.min(xy[:, 1]), np.max(xy[:, 1]))

    def add_attribute(self, name):
        if self.has_attribute(name):
            raise AttributeError(
                'Non-unique attribute name: '
                + 'Attribute attribute name already exists.')
        else:
            self.__attributes[name] = None

    def has_attribute(self, name):
        if name in self.__attributes.keys():
            return True
        else:
            return False

    def get_attribute(self, name):
        if not self.has_attribute(name):
            raise AttributeError('Attribute {} not set.'.format(name))
        return self.__attributes[name]

    def set_attribute(self, name, values, elements=False):
        if name not in self.get_attribute_names():
            raise AttributeError(
                'Cannot set attribute: {} is not an attribute.'.format(name))
        values = np.asarray(values)
        assert isinstance(elements, bool)
        if elements:
            assert values.shape[0] == self.elements.shape[0]
        else:
            assert values.shape[0] == self.vertices.shape[0]
        self.__attributes[name] = values

    def remove_attribute(self, name):
        if name in self.get_attribute_names():
            self.__attributes.pop(name)
        else:
            raise AttributeError(
                'Cannot remove attribute: attribute does not exist.')

    def get_attribute_names(self):
        return list(self.__attributes.keys())

    def interpolate(self, Dataset):
        assert isinstance(Dataset, gdal.Dataset)
        if not self.SpatialReference.IsSame(
                    gdal_tools.get_SpatialReference(Dataset)):
            Dataset = gdal_tools.Warp(Dataset, dstSRS=self.SpatialReference)
        x, y, z = gdal_tools.get_arrays(Dataset)
        bbox = gdal_tools.get_Bbox(Dataset)
        f = RectBivariateSpline(x, y, z.T, bbox=[bbox.xmin, bbox.xmax,
                                                 bbox.ymin, bbox.ymax])
        idxs = np.where(np.logical_and(
                            np.logical_and(
                                bbox.xmin <= self.vertices[:, 0],
                                bbox.xmax >= self.vertices[:, 0]),
                            np.logical_and(
                                bbox.ymin <= self.vertices[:, 1],
                                bbox.ymax >= self.vertices[:, 1])))[0]
        values = f.ev(self.vertices[idxs, 0], self.vertices[idxs, 1])
        new_values = self.values.copy()
        for i, idx in enumerate(idxs):
            new_values[idx] = values[i]
        self.values = new_values

    def has_invalid(self):
        return np.any(np.isnan(self.values))

    def fix_invalid(self, method='nearest'):
        if self.has_invalid():
            if method == 'nearest':
                idx = np.where(~np.isnan(self.values))
                _idx = np.where(np.isnan(self.values))
                values = griddata(
                    (self.x[idx], self.y[idx]), self.values[idx],
                    (self.x[_idx], self.y[_idx]), method='nearest')
                new_values = self.values.copy()
                for i, idx in enumerate(_idx):
                    new_values[idx] = values[i]
                self._values = new_values
                return self.values
            else:
                raise NotImplementedError

    def make_plot(self, show=False, levels=256):
        z = np.ma.masked_invalid(self.values)
        vmin, vmax = z.min(), z.max()
        z = z.filled(fill_value=-99999.)
        if isinstance(levels, int):
            levels = np.linspace(vmin, vmax, levels)
        plt.tricontourf(self.mpl_tri, z, levels=levels)
        plt.gca().axis('scaled')
        if show:
            plt.show()
        plt.gca().axis('scaled')
        return plt.gca()

    @staticmethod
    def parse_gr3(path):
        gr3 = dict()
        gr3['x'] = list()
        gr3['y'] = list()
        gr3['values'] = list()
        gr3['node_id'] = list()
        gr3['elements'] = list()
        gr3['element_id'] = list()
        gr3['ocean_boundaries'] = list()
        gr3['land_boundaries'] = list()
        gr3['inner_boundaries'] = list()
        gr3['inflow_boundaries'] = list()
        gr3['outflow_boundaries'] = list()
        gr3['weir_boundaries'] = list()
        gr3['culvert_boundaries'] = list()
        with open(path, 'r') as f:
            gr3['description'] = "{}".format(f.readline())
            NE, NP = map(int, f.readline().split())
            _NP = len([])
            while _NP < NP:
                node_id, x, y, z = f.readline().split()
                gr3['node_id'].append(int(node_id)-1)
                gr3['x'].append(float(x))
                gr3['y'].append(float(y))
                gr3['values'].append(-float(z))
                _NP += 1
            _NE = len([])
            while _NE < NE:
                line = f.readline().split()
                gr3['element_id'].append(int(line[0])-1)
                if int(line[1]) != 3:
                    raise NotImplementedError(
                        'Package only supports triangular meshes.')
                gr3['elements'].append([int(x)-1 for x in line[2:]])
                _NE += 1
            # Assume EOF if NOPE is empty.
            try:
                NOPE = int(f.readline().split()[0])
            except IndexError:
                return f
            # For now, let NOPE=-1 mean a self closing mesh
            # reassigning NOPE to 0 until further implementation is applied.
            if NOPE == -1:
                NOPE = 0
            _NOPE = len([])
            f.readline()  # Number of total open ocean nodes. Not used.
            while _NOPE < NOPE:
                gr3['ocean_boundaries'].append({'indexes': list()})
                NETA = int(f.readline().split()[0])
                _NETA = len([])
                while _NETA < NETA:
                    gr3['ocean_boundaries'][_NOPE]['indexes'].append(
                                                int(f.readline().split()[0])-1)
                    _NETA += 1
                _NOPE += 1
            NBOU = int(f.readline().split()[0])
            _NBOU = len([])
            f.readline()
            while _NBOU < NBOU:
                NVELL, IBTYPE = map(int, f.readline().split()[:2])
                _NVELL = 0
                if IBTYPE in [0, 10, 20]:
                    gr3['land_boundaries'].append(
                        {'ibtype': IBTYPE,
                         'indexes': list()})
                elif IBTYPE in [1, 11, 21]:
                    gr3['inner_boundaries'].append(
                        {'ibtype': IBTYPE,
                         'indexes': list()})
                elif IBTYPE in [2, 12, 22, 102, 122]:
                    gr3['inflow_boundaries'].append(
                        {'ibtype': IBTYPE,
                         'indexes': list()})
                elif IBTYPE in [3, 13, 23]:
                    gr3['outflow_boundaries'].append(
                        {'ibtype': IBTYPE,
                         'indexes': list(),
                         'barrier_heights': list(),
                         'supercritical_flow_coefficients': list()})

                elif IBTYPE in [4, 24]:
                    gr3['weir_boundaries'].append(
                        {'ibtype': IBTYPE,
                         'front_face_indexes': list(),
                         'back_face_indexes': list(),
                         'barrier_heights': list(),
                         'subcritical_flow_coefficients': list(),
                         'supercritical_flow_coefficients': list()})
                elif IBTYPE in [5, 25]:
                    gr3['culvert_boundaries'].append(
                        {'ibtype': IBTYPE,
                         'front_face_indexes': list(),
                         'back_face_indexes': list(),
                         'barrier_heights': list(),
                         'subcritical_flow_coefficients': list(),
                         'supercritical_flow_coefficients': list(),
                         'cross_barrier_pipe_heights': list(),
                         'friction_factors': list(),
                         'pipe_diameters': list()})
                else:
                    raise Exception('IBTYPE={} '.format(IBTYPE)
                                    + 'found in fort.14 not recongnized. ')
                while _NVELL < NVELL:
                    line = f.readline().split()
                    if IBTYPE in [0, 10, 20]:
                        gr3['land_boundaries'][-1][
                            'indexes'].append(int(line[0])-1)
                    elif IBTYPE in [1, 11, 21]:
                        gr3['inner_boundaries'][-1][
                            'indexes'].append(int(line[0])-1)
                    elif IBTYPE in [3, 13, 23]:
                        gr3['outflow_boundaries'][-1][
                            'indexes'].append(int(line[0])-1)
                        gr3['outflow_boundaries'][-1][
                            'external_barrier_heights'].append(float(line[1]))
                        gr3['outflow_boundaries'][-1][
                            'supercritical_flow_coefficients'].append(
                                float(line[2]))
                    elif IBTYPE in [2, 12, 22, 102, 122]:
                        gr3['iflowBoundaries'][-1][
                            'indexes'].append(int(line[0])-1)
                    elif IBTYPE in [4, 24]:
                        gr3['weir_boundaries'][-1][
                            'front_face_indexes'].append(int(line[0])-1)
                        gr3['weir_boundaries'][-1][
                            'back_face_indexes'].append(int(line[1])-1)
                        gr3['weir_boundaries'][-1][
                            'barrier_heights'].append(float(line[2]))
                        gr3['weir_boundaries'][-1][
                            'subcritical_flow_coefficients'].append(
                                float(line[3]))
                        gr3['weir_boundaries'][-1][
                            'supercritical_flow_coefficients'].append(
                                float(line[4]))
                    elif IBTYPE in [5, 25]:
                        gr3['culvert_boundaries'][-1][
                            'front_face_indexes'].append(int(line[0])-1)
                        gr3['culvert_boundaries'][-1][
                            'back_face_indexes'].append(int(line[1])-1)
                        gr3['culvert_boundaries'][-1][
                            'barrier_heights'].append(float(line[2]))
                        gr3['culvert_boundaries'][-1][
                            'subcritical_flow_coefficients'].append(
                                float(line[3]))
                        gr3['culvert_boundaries'][-1][
                            'supercritical_flow_coefficients'].append(
                                float(line[4]))
                        gr3['culvert_boundaries'][-1][
                            'friction_factors'].append(float(line[5]))
                        gr3['culvert_boundaries'][-1][
                            'pipe_diameters'].append(float(line[6]))
                    else:
                        Exception("Duck-typing error. "
                                  + "This exception should be unreachable.")
                    _NVELL += 1
                _NBOU += 1
        return gr3

    @property
    def vertices(self):
        return self.__vertices

    @property
    def elements(self):
        return self.__elements

    @property
    def values(self):
        return self.__values

    @property
    def x(self):
        return self.vertices[:, 0]

    @property
    def y(self):
        return self.vertices[:, 1]

    @property
    def xy(self):
        return self.vertices

    @property
    def xyz(self):
        return self.get_xyz()

    @property
    def SpatialReference(self):
        return self._get_spatial_reference()

    @property
    def planar_straight_line_graph(self):
        """
        Slow alogirthm, open for suggestions...
        """
        try:
            return self.__planar_straight_line_graph
        except AttributeError:
            _pslg = _PlanarStraightLineGraph(self.SpatialReference)
            idxs = np.vstack(list(np.where(self.mpl_tri.neighbors == -1))).T
            unique_edges = list()
            for i, j in idxs:
                unique_edges.append((self.mpl_tri.triangles[i, j],
                                     self.mpl_tri.triangles[i, (j+1) % 3]))
            unique_edges = np.asarray(unique_edges)
            ring_collection = list()
            initial_idx = 0
            for i in range(1, len(unique_edges)-1):
                if unique_edges[i-1, 1] != unique_edges[i, 0]:
                    try:
                        idx = np.where(
                            unique_edges[i-1, 1] == unique_edges[i:, 0])[0][0]
                        unique_edges[[i, idx+i]] = unique_edges[[idx+i, i]]
                    except IndexError:
                        ring_collection.append(unique_edges[initial_idx:i, 0])
                        initial_idx = i
                        continue
            if len(ring_collection) == 0:
                ring_collection.append(unique_edges[initial_idx:i, 0])
            print(
                [Path(np.hstack([self.x[ring], self.y[ring]]).T, closed=True)
                 for ring in ring_collection])
            BREAKME


            # linear_ring_collection = list()
            # for ring in ring_collection:
            #     # plt.plot()
            #     _geom = ogr.Geometry(ogr.wkbLinearRing)
            #     _geom.AssignSpatialReference(self.SpatialReference)
            #     for idx in ring:
            #         _geom.AddPoint(self.x[idx], self.y[idx])
            #     _geom.CloseRings()
            #     _polygon = ogr.Geometry(ogr.wkbPolygon)
            #     _polygon.AddGeometry(_geom)
            #     linear_ring_collection.append(_polygon)
            # idx = np.where(
            #     np.max([geom.GetArea() for geom in linear_ring_collection]))[0]
            # linear_ring = linear_ring_collection.pop(idx[0])
            # # contains = [linear_ring.Contains(geom)
            # #             for geom in linear_ring_collection]
            # contains = [geom.Within(linear_ring)
            #             for geom in linear_ring_collection]
            # print(contains)
            # BREAKME

                # print(_polygon.IsValid())
                # _geom.CloseRings()
                # print(_geom.IsValid())
                # BREAKME
                # vertices = np.asarray([(x, y) for x, y, _ in _geom.GetPoints()])
                # plt.plot(vertices[:, 0], vertices[:, 1])
                # plt.show()

                # _geom.AssignSpatialReference(self.SpatialReference)
                # print(_geom.IsValid())
                # BREAKME
                # linear_ring_collection.append(_geom)
            # print([geom.Polygonize() for geom in linear_ring_collection])
            # print([geom.IsRing() for geom in linear_ring_collection])
            # BREAKME
            # idx = np.where(
            #     np.max([geom.GetArea() for geom in linear_ring_collection]))[0]
            # linear_ring = linear_ring_collection.pop(idx[0])
            # contains = [linear_ring.Contains(geom)
            #             for geom in linear_ring_collection]
            # print(contains)
            # BREAKME

        # print(linear_ring)

        # while len(linear_ring_collection) > 0:
        #     linear_ring = linear_ring_collection.pop()
        #     for remaining_rings in linear_ring_collection:


                



        # print(geom_collection)
        # BREAKME

        # ring_collection = []
        # while len(ring_collection) > 0:
        #     ring = ring_collection.pop()
        #     for test_ring in ring_collection:
        #         if ring.contains(test_ring)
        # multipolygon = ogr.Geometry(ogr.wkbMultiPolygon)
        # for polygon in polygons:
        #     multipolygon.AddGeometry(polygon)
        # #  -------------------
        # geom_collection = list()
        # for ring in ring_collection:
        #     _geom = ogr.Geometry(ogr.wkbLinearRing)
        #     _geom.AssignSpatialReference(self.SpatialReference)
        #     for idx in ring:
        #         _geom.AddPoint_2D(self.x[idx[0]], self.y[idx[0]])
        #     geom_collection.append(_geom)
        # lengths = [_geom.Length() for _geom in geom_collection]
        # outer_edges = ring_collection.pop(
        #         np.where(np.max(lengths) == lengths)[0][0])
        # inner_edges = ring_collection
        # outer_vertices = self.vertices[outer_edges[:, 0]]
        # outer_vertices = np.vstack([outer_vertices, outer_vertices[0, :]])
        # inner_vertices = [self.vertices[ring[:, 0]] for ring in inner_edges]
        # inner_vertices = [np.vstack([vertices, vertices[0, :]])
        #                   for vertices in inner_vertices]
        # return _PlanarStraightLineGraph(
        #         self.SpatialReference, outer_vertices, *inner_vertices,
        #         outer_edges=outer_edges, inner_edges=inner_edges)
            return self.__planar_straight_line_graph

    # def make_plot(self):
    #     plt.plot(self.outer_vertices[:, 0], self.outer_vertices[:, 1])
    #     for vertices in self.inner_vertices:
    #         plt.plot(vertices[:, 0], vertices[:, 1])
    #     plt.show()

    @property
    def mpl_tri(self):
        # Since vertices, elements are constants, Triangulation cannot change.
        try:
            return self.__mpl_tri
        except AttributeError:
            if not np.ma.is_masked(self.values):
                self.__mpl_tri = Triangulation(self.x, self.y, self.elements)
                return self.__mpl_tri
            else:
                self.__mpl_tri = Triangulation(
                    self.x, self.y, self.elements,
                    mask=np.any(self.values.mask[self.elements], axis=1))
                return self.__mpl_tri

    @property
    def ndim(self):
        return 2

    @property
    def num_elements(self):
        return self.elements.shape[0]

    @property
    def num_nodes(self):
        return self.vertices.shape[0]

    @SpatialReference.setter
    def SpatialReference(self, SpatialReference):
        assert isinstance(SpatialReference, (int, osr.SpatialReference)), \
            "Input must be a EPSG code or osr.SpatialReference instance."
        msg = "Mesh must have a spatial reference assigned before "
        msg += "transformation can occur."
        assert self._get_spatial_reference() is not None, msg
        self._vertices = self.transform_vertices(
                self.__vertices, self._get_spatial_reference(),
                SpatialReference)
        self._SpatialReference = SpatialReference

    @property
    def _vertices(self):
        return self.__vertices

    @property
    def _elements(self):
        return self.__elements

    @property
    def _values(self):
        return self.__values

    @property
    def _SpatialReference(self):
        return self._get_spatial_reference()

    @_vertices.setter
    def _vertices(self, vertices):
        vertices = np.asarray(vertices)
        assert vertices.shape[1] == self.ndim
        self.__vertices = vertices

    @_elements.setter
    def _elements(self, elements):
        elements = np.asarray(elements)
        assert elements.shape[1] == 3
        self.__elements = elements

    @_values.setter
    def _values(self, values):
        if values is None:
            values = np.full((self.vertices.shape[0],), np.nan)
        if not np.ma.is_masked(values):
            values = np.asarray(values)
        assert values.shape[0] == self.vertices.shape[0]
        self.__values = values

    @_SpatialReference.setter
    def _SpatialReference(self, SpatialReference):
        if SpatialReference is not None:
            self._set_spatial_reference(SpatialReference)
        else:
            self._clear_spatial_reference()
