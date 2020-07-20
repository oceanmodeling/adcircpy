from collections import defaultdict
from collections.abc import Mapping
from functools import lru_cache
from itertools import permutations
import logging
import pathlib
import uuid

from haversine import Unit, haversine
from matplotlib.collections import PolyCollection
from matplotlib.transforms import Bbox
from matplotlib.tri import Triangulation
import numpy
import numpy as np
from pyproj import CRS, Proj, Transformer
from shapely.geometry import MultiLineString, MultiPolygon, Polygon
from shapely.ops import polygonize, unary_union

from adcircpy.mesh import grd, sms2dm
from adcircpy.mesh.figures import _figure as _fig


class EuclideanMesh2D:
    def __init__(
        self,
        coords,
        triangles=None,
        quads=None,
        values=None,
        crs=None,
        description=None,
    ):
        self._coords = coords
        self._triangles = triangles
        self._quads = quads
        self._values = values
        self._crs = crs
        self._description = description

    @classmethod
    def open(cls, path, crs=None, fmt="grd"):
        assert fmt.lower() in [
            'grd', 'gr3', 'adcirc', 'schism',
            # '2dm', 'sms',
            # 'msh'
        ]

        if fmt.lower() in ['grd', 'gr3', 'adcirc', 'schism']:
            return cls.open_grd(path, crs)

        # elif fmt.lower() in ['2dm', 'sms']:
        #     return cls.open_2dm(path, crs)

    # @classmethod
    # def open_2dm(cls, path, crs=None):
    #     raise NotImplementedError

    @classmethod
    def open_grd(cls, path, crs=None):
        grid = grd.reader(path)
        grid.update({"crs": crs})
        return cls.from_grd(grid)

    @classmethod
    def open_gr3(cls, path, crs=None):
        return cls.open_grd(path, crs)

    @classmethod
    def from_grd(cls, grid):
        """
        grd is a dictionary of of the form:
        """
        return cls(**grd.euclidean_mesh(grid))

    def transform_to(self, dst_crs):
        dst_crs = CRS.from_user_input(dst_crs)
        if self.srs != dst_crs.srs:
            transformer = Transformer.from_crs(
                self.crs, dst_crs,
                always_xy=True
            )
            xy = list(zip(*transformer.transform(self.x, self.y)))
            ids = list(self._coords.keys())
            self._coords = {ids[i]: coord for i, coord in enumerate(xy)}
            self._crs = dst_crs

    def get_node_index(self, id):
        return self.node_index[id]

    def get_node_id(self, index):
        return self.node_id[index]

    def get_element_index(self, id):
        return self.element_index[id]

    def get_element_id(self, index):
        return self.element_id[index]

    def get_x(self, crs=None):
        return self.get_xy(crs)[:, 0]

    def get_y(self, crs=None):
        return self.get_xy(crs)[:, 1]

    def get_xy(self, crs=None):
        if crs is not None:
            crs = CRS.from_user_input(crs)
            if crs.srs != self.srs:
                transformer = Transformer.from_crs(
                    self.crs, crs, always_xy=True)
                x, y = transformer.transform(self.x, self.y)
                return np.vstack([x, y]).T
        return np.vstack([self.x, self.y]).T

    def get_extent(self, crs=None):
        xy = self.get_xy(crs)
        return (np.min(xy[:, 0]), np.max(xy[:, 0]),
                np.min(xy[:, 1]), np.max(xy[:, 1]))

    def get_bbox(self, crs=None):
        xmin, xmax, ymin, ymax = self.get_extent(crs)
        return Bbox([[xmin, ymin], [xmax, ymax]])

    def write(self, path, overwrite=False, fmt='gr3'):
        path = pathlib.Path(path)
        if path.is_file() and not overwrite:
            msg = 'File exists, use overwrite=True to allow overwrite.'
            raise IOError(msg)
        with open(path, 'w') as f:
            f.write(self.ascii_string(fmt))

    def ascii_string(self, fmt):
        if fmt.lower() in ['grd', 'gr3', 'adcirc', 'schism']:
            return self.grd
        elif fmt.lower() in ['sms', '2dm']:
            return self.sms2dm
        else:
            msg = f"File format {fmt} not recognized."
            raise Exception(msg)

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
        msg = "values cannot be None"
        assert values is not None, msg
        values = np.array(values)
        assert isinstance(elements, bool)
        if elements:
            assert values.shape[0] == self.elements.shape[0]
        else:
            assert values.shape[0] == self.coords.shape[0]
        self._attributes[name]['values'] = values
        self._attributes[name].update(properties)

    def remove_attribute(self, name):
        if name in self.get_attribute_names():
            self._attributes.pop(name)
        else:
            raise AttributeError(
                'Cannot remove attribute: attribute does not exist.')

    @_fig
    def tricontourf(self, axes=None, show=True, figsize=None, **kwargs):
        if len(self.tria3) > 0:
            axes.tricontourf(self.triangulation, self.values, **kwargs)
        return axes

    @_fig
    def tripcolor(self, axes=None, show=True, figsize=None, **kwargs):
        if len(self.tria3) > 0:
            axes.tripcolor(self.triangulation, self.values, **kwargs)
        return axes

    @_fig
    def triplot(
        self,
        axes=None,
        show=False,
        figsize=None,
        linewidth=0.07,
        color='black',
        **kwargs
    ):
        if len(self.triangles) > 0:
            kwargs.update({'linewidth': linewidth})
            kwargs.update({'color': color})
            axes.triplot(self.triangulation, **kwargs)
        return axes

    @_fig
    def quadplot(
        self,
        axes=None,
        show=False,
        figsize=None,
        facecolor='none',
        edgecolor='k',
        linewidth=0.07,
        **kwargs
    ):
        if len(self.quads) > 0:
            pc = PolyCollection(
                self.coords[self.quads],
                facecolor=facecolor,
                edgecolor=edgecolor,
                linewidth=0.07,
            )
            axes.add_collection(pc)
        return axes

    @_fig
    def quadface(
        self,
        axes=None,
        show=False,
        figsize=None,
        **kwargs
    ):
        if len(self.quad4) > 0:
            pc = PolyCollection(
                self.coords[self.quad4],
                **kwargs
            )
            quad_value = np.mean(self.values[self.quad4], axis=1)
            pc.set_array(quad_value)
            axes.add_collection(pc)
        return axes

    @_fig
    def plot_wireframe(self, axes=None, show=False, **kwargs):
        axes = self.triplot(axes=axes, **kwargs)
        axes = self.quadplot(axes=axes, **kwargs)
        return axes

    @property
    @lru_cache(maxsize=None)
    def coords(self):
        return np.array(
            [coord for coord in self._coords.values()]
        )

    @property
    def xy(self):
        return self.coords

    @property
    def x(self):
        return self.coords[:, 0]

    @property
    def y(self):
        return self.coords[:, 1]

    @property
    @lru_cache(maxsize=None)
    def values(self):
        return self._values

    @property
    def triangulation(self):
        if len(self.tria3) > 0:
            return Triangulation(self.x, self.y, self.tria3)

    @property
    def triangles(self):
        return self.tria3

    @property
    def quads(self):
        return self.quad4

    @property
    @lru_cache(maxsize=None)
    def nodes(self):
        return list(self._nodes.items())

    @property
    @lru_cache(maxsize=None)
    def elements(self):
        return list(self._elements.values())

    @property
    @lru_cache(maxsize=None)
    def coords_id(self):
        return self._coords.keys()

    @property
    @lru_cache(maxsize=None)
    def triangles_id(self):
        return self._triangles.keys()

    @property
    @lru_cache(maxsize=None)
    def quads_id(self):
        return self._quads.keys()

    @property
    @lru_cache(maxsize=None)
    def node_id(self):
        return {index: id for index, id in enumerate(self._coords)}

    @property
    @lru_cache(maxsize=None)
    def element_id(self):
        return {index: id for index, id in enumerate(self._elements)}

    @property
    @lru_cache(maxsize=None)
    def node_index(self):
        return {id: index for index, id in enumerate(self._coords)}

    @property
    @lru_cache(maxsize=None)
    def element_index(self):
        return {id: index for index, id in enumerate(self._elements)}

    @property
    def bbox(self):
        return self.get_bbox()

    @property
    def description(self):
        return self._description

    @property
    def proj(self):
        return Proj(self.crs)

    @property
    def srs(self):
        return self.proj.srs

    @property
    def crs(self):
        return self._crs

    @property
    def grd(self):
        return grd.string(self._grd)

    @property
    def sms2dm(self):
        return sms2dm.string(self._sms2dm)

    @property
    @lru_cache(maxsize=None)
    def tria3(self):
        return np.array(
            [list(map(self.get_node_index, index))
             for index in self._triangles.values()])

    @property
    @lru_cache(maxsize=None)
    def quad4(self):
        return np.array(
            [list(map(self.get_node_index, index))
             for index in self._quads.values()])

    @lru_cache(maxsize=None)
    def mesh_hull(self, element_lengths=None):
        if element_lengths is None:
            elements = list(self.triangles) + list(self.quads)
        else:
            elements = []
            if 3 in element_lengths:
                elements.extend(self.triangles)
            elif 4 in element_lengths:
                elements.extend(self.quads)
        return polygon_from_vertices(boundary_points(self.coords, elements))

    @property
    @lru_cache(maxsize=None)
    def index_ring_collection(self):
        # find boundary edges using triangulation neighbors table,
        # see: https://stackoverflow.com/a/23073229/7432462

        # TODO add in quads (by setting `element_lengths` to `None`)
        mesh_hull = self.mesh_hull(element_lengths=(3,))

        index_ring_collection = {
            index: {
                'exterior' : polygon.exterior,
                'interiors': [interior for interior in polygon.interiors]
            } for index, polygon in enumerate(mesh_hull)
        }

        for index, ring in index_ring_collection.items():
            exterior = index_ring_collection[index]['exterior'].coords
            interiors = [interior.coords for interior in
                         index_ring_collection[index]['interiors']]

            exterior = [(exterior[index],
                         exterior[index - len(exterior) + 1])
                        for index in range(len(exterior))]
            interiors = [[(interior[index],
                           interior[index - len(interior) + 1])
                          for index in range(len(interior))]
                         for interior in interiors]

            index_ring_collection[index]['exterior'] = numpy.array([[
                [tuple(coord) for coord in self.coords].index(vertex)
                for vertex in edge] for edge in exterior])
            index_ring_collection[index]['interiors'] = [
                numpy.array([[
                    [tuple(coord) for coord in self.coords].index(vertex)
                    for vertex in edge] for edge in interior])
                for interior in interiors
            ]

        return index_ring_collection

    @property
    @lru_cache(maxsize=None)
    def outer_ring_collection(self):
        outer_ring_collection = defaultdict()
        for key, ring in self.index_ring_collection.items():
            outer_ring_collection[key] = ring['exterior']
        return outer_ring_collection

    @property
    @lru_cache(maxsize=None)
    def inner_ring_collection(self):
        inner_ring_collection = defaultdict()
        for key, rings in self.index_ring_collection.items():
            inner_ring_collection[key] = rings['interiors']
        return inner_ring_collection

    @property
    @lru_cache(maxsize=None)
    def node_neighbors(self):
        node_neighbors = defaultdict(set)
        for simplex in self.triangulation.triangles:
            for i, j in permutations(simplex, 2):
                node_neighbors[i].add(j)
        return node_neighbors

    @property
    @lru_cache(maxsize=None)
    def node_distances_meters(self):
        points = self.get_xy("EPSG:4326")
        node_distances = {}
        for k, v in self.node_neighbors.items():
            x0, y0 = points[k]
            node_distances[k] = {}
            for idx in v:
                x1, y1 = points[idx]
                node_distances[k][idx] = haversine(
                    (x0, y0),
                    (x1, y1),
                    unit=Unit.METERS
                )
        return node_distances

    @property
    @lru_cache(maxsize=None)
    def _logger(self):
        return logging.getLogger(__name__ + '.' + self.__class__.__name__)

    @description.setter
    def description(self, description):
        self._description = description

    def _certify_input_geom(self, geom_type, geom):
        geom_types = {
            "triangles": 3,
            "quads"    : 4,
        }
        _geom = np.asarray(list(geom.values()))
        if len(_geom) > 0:
            msg = f'Invalid shape for geom {_geom.shape}.'
            assert _geom.shape[1] == geom_types[geom_type], msg
            msg = f"{geom_type} members must be a subset of coords keys"
            for IDtags in geom.values():
                for IDtag in IDtags:
                    assert IDtag in self.coords_id, msg

    @property
    def _coords(self):
        return self.__coords

    @_coords.setter
    def _coords(self, coords):
        msg = "coord argument must be a dictionary of the form "
        msg += "\\{coord_id:  (x, y)\\}"
        assert isinstance(coords, Mapping), msg
        for coord in coords.values():
            assert len(coord) == 2, msg
            assert isinstance(coord[0], (float, int)), msg
            assert isinstance(coord[1], (float, int)), msg
        self.__coords = coords
        type(self).coords.fget.cache_clear()

    @property
    def _triangles(self):
        return self.__triangles

    @property
    def _quads(self):
        return self.__quads

    @property
    def _values(self):
        return self.__values

    @property
    def _crs(self):
        return self.__crs

    @property
    def _description(self):
        return self.__description

    @property
    @lru_cache(maxsize=None)
    def _nodes(self):
        return {id: ((x, y), -self.values[i]) for i, (id, (x, y))
                in enumerate(self._coords.items())}

    @property
    @lru_cache(maxsize=None)
    def _elements(self):
        elements_id = list()
        elements_id.extend(list(self._triangles.keys()))
        elements_id.extend(list(self._quads.keys()))
        elements_id = range(1, len(elements_id) + 1) \
            if len(set(elements_id)) != len(elements_id) else elements_id
        elements = list()
        elements.extend(list(self._triangles.values()))
        elements.extend(list(self._quads.values()))
        elements = {
            elements_id[i]: indexes for i, indexes in enumerate(elements)}
        return elements

    @property
    @lru_cache(maxsize=None)
    def _grd(self):
        description = self.description
        if self.crs is not None and self.crs.srs not in description:
            description += f"; {self.crs.srs}"
        return {
            "nodes"      : self._nodes,
            "elements"   : self._elements,
            "description": description,
        }

    @property
    @lru_cache(maxsize=None)
    def _sms2dm(self):
        description = self.description
        if self.crs is not None:
            description += f"; {self.crs.srs}"
        return {
            "ND" : self._nodes,
            "E3T": self._triangles,
            "E4Q": self._quads,
        }

    @_values.setter
    def _values(self, values):
        if values is None:
            values = []
        values = np.asarray(values)
        if len(values) > 0:
            if len(values.shape) in [1, 2]:
                assert values.shape[0] == self.coords.shape[0]
            elif len(values.shape) == 3:
                assert values.shape[1] == self.coords.shape[0]
            else:
                msg = f"input values has invalid shape: {values.shape}"
                raise Exception(msg)

        else:
            values = np.full(self.coords.shape[0], np.nan)
        self.__values = values
        # type(self).values.fget.cache_clear()

    @_triangles.setter
    def _triangles(self, triangles):
        if triangles is None:
            triangles = {}
        self._certify_input_geom("triangles", triangles)
        self.__triangles = triangles

    @_quads.setter
    def _quads(self, quads):
        if quads is None:
            quads = {}
        self._certify_input_geom("quads", quads)
        self.__quads = quads

    @_crs.setter
    def _crs(self, crs):
        if crs is not None:
            crs = CRS.from_user_input(crs)
        self.__crs = crs

    @_description.setter
    def _description(self, description):
        if description is None:
            description = uuid.uuid4().hex[:8]
        assert isinstance(description, str)
        self.__description = description

    # auxilliary functions
    @staticmethod
    def sort_edges(edges):
        if len(edges) == 0:
            return edges
        # start ordering the edges into linestrings
        edge_collection = list()
        ordered_edges = [edges.pop(-1)]
        e0, e1 = [list(t) for t in zip(*edges)]
        while len(edges) > 0:
            if ordered_edges[-1][1] in e0:
                idx = e0.index(ordered_edges[-1][1])
                ordered_edges.append(edges.pop(idx))
            elif ordered_edges[0][0] in e1:
                idx = e1.index(ordered_edges[0][0])
                ordered_edges.insert(0, edges.pop(idx))
            elif ordered_edges[-1][1] in e1:
                idx = e1.index(ordered_edges[-1][1])
                ordered_edges.append(
                    list(reversed(edges.pop(idx))))
            elif ordered_edges[0][0] in e0:
                idx = e0.index(ordered_edges[0][0])
                ordered_edges.insert(
                    0, list(reversed(edges.pop(idx))))
            else:
                edge_collection.append(tuple(ordered_edges))
                idx = -1
                ordered_edges = [edges.pop(idx)]
            e0.pop(idx)
            e1.pop(idx)
        # finalize
        if len(edge_collection) == 0 and len(edges) == 0:
            edge_collection.append(tuple(ordered_edges))
        else:
            edge_collection.append(tuple(ordered_edges))
        return edge_collection

    @staticmethod
    def signed_polygon_area(vertices):
        # https://code.activestate.com/recipes/578047-area-of-polygon-using-shoelace-formula/
        n = len(vertices)  # of vertices
        area = 0.0
        for i in range(n):
            j = (i + 1) % n
            area += vertices[i][0] * vertices[j][1]
            area -= vertices[j][0] * vertices[i][1]
        return area / 2.0

    @property
    @lru_cache(maxsize=None)
    def _attributes(self):
        return {}


def boundary_indices(elements: [[int]]) -> numpy.array:
    shape_lengths = {len(element) for element in elements}

    elements = [[element for element in elements
                 if len(element) == shape_length]
                for shape_length in shape_lengths]

    boundary_edge_indices = []
    for shape_elements in elements:
        shape_elements = numpy.stack(shape_elements, axis=0)
        shape_length = shape_elements.shape[1]
        indices = numpy.stack([(shape_elements[:, index],
                                shape_elements[:, index - shape_length + 1])
                               for index in range(shape_length)], axis=1).T

        indices = numpy.reshape(
            indices, (indices.shape[0] * indices.shape[1], indices.shape[2])
        )

        unique_indices, counts = numpy.unique(indices, axis=0,
                                              return_counts=True)
        boundary_edge_indices.extend(unique_indices[counts == 1])

    return boundary_edge_indices


def boundary_points(points: numpy.array, elements: [[int]]) -> numpy.array:
    return numpy.stack(points[boundary_indices(elements)], axis=0)


def polygon_from_vertices(vertices: numpy.array) -> MultiPolygon:
    polygons = polygonize(MultiLineString(vertices.tolist()))
    polygons = remove_interior_duplicates(polygons)
    return unary_union(polygons)


def remove_interior_duplicates(polygons: [Polygon]) -> [Polygon]:
    """
    Given a list of polygons, return a list without polygons whose exteriors represent the interior of another polygon.

    Parameters
    ----------
    polygons
        list of polygons

    Returns
    -------
    [Polygon]
        polygons without duplicates of interiors
    """

    if type(polygons) is not list:
        polygons = list(polygons)

    for greater_polygon in polygons:
        for interior in greater_polygon.interiors:
            interior_polygon = Polygon(interior)
            duplicate_polygons = []

            for polygon_index in range(len(polygons)):
                polygon = polygons[polygon_index]
                if polygon.equals(interior_polygon):
                    duplicate_polygons.append(polygon)

            for polygon in duplicate_polygons:
                polygons.remove(polygon)

    return polygons
