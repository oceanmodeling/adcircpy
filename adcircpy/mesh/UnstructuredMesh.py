import numpy as np
from osgeo import osr, gdal, ogr
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
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
        xy = self.transform_vertices(self.xy, self.SpatialReference,
                                     SpatialReference)
        return np.hstack([xy, self.values])

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
        self._values = new_values

    def get_planar_straight_line_graph(self):
        mpl_tri = self.mpl_tri
        idxs = np.vstack(list(np.where(mpl_tri.neighbors == -1))).T
        unique_edges = list()
        for i, j in idxs:
            unique_edges.append((mpl_tri.triangles[i, j],
                                 mpl_tri.triangles[i, (j+1) % 3]))
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
                    ring_collection.append(unique_edges[initial_idx:i, :])
                    initial_idx = i
                    continue
        if len(ring_collection) == 0:
            ring_collection.append(np.asarray(unique_edges))
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

    @property
    def vertices(self):
        return self._vertices

    @property
    def elements(self):
        return self._elements

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
    def mpl_tri(self):
        if not hasattr(self, "__mpl_tri"):
            self.__mpl_tri = Triangulation(self.x, self.y, self.elements)
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

    @values.setter
    def values(self, values):
        if self.has_attribute("__values"):
            self.remove_attribute("__values")
        self._values = values

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
        values = np.asarray(values)
        assert values.shape[0] == self.vertices.shape[0]
        self.__values = values

    @_SpatialReference.setter
    def _SpatialReference(self, SpatialReference):
        if SpatialReference is not None:
            self._set_spatial_reference(SpatialReference)
        else:
            self._clear_spatial_reference()
