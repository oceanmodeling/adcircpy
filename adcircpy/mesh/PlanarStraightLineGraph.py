import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.path import Path as _Path
from osgeo import osr, ogr
# from adcircpy.lib._SpatialReference import _SpatialReference


class PlanarStraightLineGraph:
    """
    Planar straight line graph representing the outer and inner hulls of a
    geospatial domain.
    """

    def __init__(self, SpatialReference):
        self.__geom = ogr.Geometry(ogr.wkbMultiPolygon)

    def add_polygon(
        self,
        outer_vertices,
        *inner_vertices,
        SpatialReference=None
    ):
        vertices = np.asarray(vertices)
        for x, y in vertices:
            self.__geom.AddPoint_2D(x, y)

    def add_ogrPolygon(self, ogrPolygon):
        pass

    @classmethod
    def from_ogrMultiPolygon(cls, ogrMultiPolygon):
        raise NotImplementedError

    @property
    def SpatialReference(self):
        return self._SpatialReference

    @SpatialReference.setter
    def SpatialReference(self, SpatialReference):
        if isinstance(SpatialReference, int):
            EPSG = SpatialReference
            SpatialReference = osr.SpatialReference()
            SpatialReference.ImportFromEPSG(EPSG)
        assert isinstance(SpatialReference, osr.SpatialReference)
        if not self.SpatialReference.IsSame(SpatialReference):
            CoordinateTransform = osr.CoordinateTransformation(
                                                        self.SpatialReference,
                                                        SpatialReference)
            # convert outer vertices
            if self.SpatialReference.IsGeographic():
                vertices = [(y, x) for x, y in self.outer_vertices]
            else:
                vertices = [(x, y) for x, y in self.outer_vertices]
            vertices = CoordinateTransform.TransformPoints(vertices)
            vertices = np.asarray([(x, y) for x, y, _ in vertices])
            self._outer_vertices = vertices
            inner_vertices = list()
            for vertices in self.inner_vertices:
                if self.SpatialReference.IsGeographic():
                    vertices = [(y, x) for x, y in vertices]
                else:
                    vertices = [(x, y) for x, y in vertices]
                vertices = CoordinateTransform.TransformPoints(vertices)
                inner_vertices.append(
                    np.asarray([(x, y) for x, y, _ in vertices]))
            self._inner_vertices = inner_vertices
            self._SpatialReference = SpatialReference




    # @property
    # def outer_vertices(self):
    #     return self._outer_vertices

    # @property
    # def inner_vertices(self):
    #     return self._inner_vertices

    # @property
    # def outer_edges(self):
    #     return self._outer_edges

    # @property
    # def inner_edges(self):
    #     return self._inner_edges

    # @property
    # def ndim(self):
    #     return 2

    # @property
    # def vert2(self):
    #     vert2 = self.outer_vertices
    #     for i, array in enumerate(self.inner_vertices):
    #         vert2 = np.vstack([vert2, array])
    #     return vert2

    # @property
    # def edge2(self):
    #     edge2 = self.outer_wire_network_edges
    #     i = edge2[-1, 0] + 1
    #     for _i, _ in enumerate(self.inner_vertices):
    #         edge2 = np.vstack([edge2, self.inner_wire_network_edges[_i]+i])
    #         i = edge2[-1, 0] + 1
    #     return edge2

    # @property
    # def outer_wire_network_edges(self):
    #     if not hasattr(self, "__outer_wire_network_edges"):
    #         outer_edges = list()
    #         for i, _ in enumerate(self.outer_vertices):
    #             if i != len(self.outer_vertices)-1:
    #                 outer_edges.append((i, i+1))
    #             else:
    #                 outer_edges.append((i, 0))
    #         self.__outer_wire_network_edges = np.asarray(outer_edges)
    #     return self.__outer_wire_network_edges

    # @property
    # def inner_wire_network_edges(self):
    #     if not hasattr(self, "__inner_wire_network_edges"):
    #         inner_edges = list()
    #         for inner_vertices in self.inner_vertices:
    #             _inner_edges = list()
    #             for i, _ in enumerate(inner_vertices):
    #                 if i != len(inner_vertices)-1:
    #                     _inner_edges.append((i, i+1))
    #                 else:
    #                     _inner_edges.append((i, 0))
    #             inner_edges.append(np.asarray(_inner_edges))
    #         self.__inner_wire_network_edges = inner_edges
    #     return self.__inner_wire_network_edges

    # @property
    # def _outer_vertices(self):
    #     return self.__outer_vertices

    # @property
    # def _inner_vertices(self):
    #     return self.__inner_vertices

    # @property
    # def _outer_edges(self):
    #     return self.__outer_edges

    # @property
    # def _inner_edges(self):
    #     return self.__inner_edges

    # @_outer_vertices.setter
    # def _outer_vertices(self, outer_vertices):
    #     outer_vertices = np.asarray(outer_vertices)
    #     assert outer_vertices.shape[1] == 2
    #     self.__outer_vertices = outer_vertices

    # @_inner_vertices.setter
    # def _inner_vertices(self, inner_vertices):
    #     inner_vertices = [np.asarray(_) for _ in inner_vertices]
    #     if len(inner_vertices) > 0:
    #         outer_hull = _Path(self.outer_vertices, closed=True)
    #         for coords in inner_vertices:
    #             assert np.all(outer_hull.contains_points(coords))
    #     self.__inner_vertices = inner_vertices

    # @_outer_edges.setter
    # def _outer_edges(self, outer_edges):
    #     if outer_edges is not None:
    #         outer_edges = np.asarray(outer_edges)
    #         assert outer_edges.shape[1] == self.outer_vertices.shape[1]
    #     else:
    #         outer_edges = self.outer_wire_network_edges
    #     self.__outer_edges = outer_edges

    # @_inner_edges.setter
    # def _inner_edges(self, inner_edges):
    #     if inner_edges is not None:
    #         assert len(inner_edges) == len(self.inner_vertices)
    #         inner_edges = [np.asarray(_) for _ in inner_edges]
    #         for i, array in enumerate(inner_edges):
    #             assert array.shape[1] == self.inner_vertices[i].shape[1]
    #     else:
    #         inner_edges = self.inner_wire_network_edges
    #     self.__inner_edges = inner_edges
