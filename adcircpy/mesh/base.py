import logging
import uuid
import pathlib
from functools import lru_cache
from collections.abc import Mapping
from collections import defaultdict
import numpy as np
from matplotlib.tri import Triangulation
from matplotlib.path import Path
from matplotlib.collections import PolyCollection
from matplotlib.transforms import Bbox
from shapely.geometry import Polygon
from pyproj import Proj, CRS, Transformer
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

    @property
    @lru_cache(maxsize=None)
    def index_ring_collection(self):
        # find boundary edges using triangulation neighbors table,
        # see: https://stackoverflow.com/a/23073229/7432462
        boundary_edges = list()
        tri = self.triangulation
        idxs = np.vstack(
            list(np.where(tri.neighbors == -1))).T
        for i, j in idxs:
            boundary_edges.append(
                (int(tri.triangles[i, j]),
                    int(tri.triangles[i, (j+1) % 3])))
        index_ring_collection = self.sort_edges(boundary_edges)
        # sort index_rings into corresponding "polygons"
        areas = list()
        vertices = self.coords
        for index_ring in index_ring_collection:
            e0, e1 = [list(t) for t in zip(*index_ring)]
            areas.append(float(Polygon(vertices[e0, :]).area))

        # maximum area must be main mesh
        idx = areas.index(np.max(areas))
        exterior = index_ring_collection.pop(idx)
        areas.pop(idx)
        _id = 0
        _index_ring_collection = dict()
        _index_ring_collection[_id] = {
            'exterior': np.asarray(exterior),
            'interiors': []
            }
        e0, e1 = [list(t) for t in zip(*exterior)]
        path = Path(vertices[e0 + [e0[0]], :], closed=True)
        while len(index_ring_collection) > 0:
            # find all internal rings
            potential_interiors = list()
            for i, index_ring in enumerate(index_ring_collection):
                e0, e1 = [list(t) for t in zip(*index_ring)]
                if path.contains_point(vertices[e0[0], :]):
                    potential_interiors.append(i)
            # filter out nested rings
            real_interiors = list()
            for i, p_interior in reversed(list(enumerate(potential_interiors))):
                _p_interior = index_ring_collection[p_interior]
                check = [index_ring_collection[_]
                         for j, _ in reversed(list(enumerate(potential_interiors)))
                         if i != j]
                has_parent = False
                for _path in check:
                    e0, e1 = [list(t) for t in zip(*_path)]
                    _path = Path(vertices[e0 + [e0[0]], :], closed=True)
                    if _path.contains_point(vertices[_p_interior[0][0], :]):
                        has_parent = True
                if not has_parent:
                    real_interiors.append(p_interior)
            # pop real rings from collection
            for i in reversed(sorted(real_interiors)):
                _index_ring_collection[_id]['interiors'].append(
                    np.asarray(index_ring_collection.pop(i)))
                areas.pop(i)
            # if no internal rings found, initialize next polygon
            if len(index_ring_collection) > 0:
                idx = areas.index(np.max(areas))
                exterior = index_ring_collection.pop(idx)
                areas.pop(idx)
                _id += 1
                _index_ring_collection[_id] = {
                    'exterior': np.asarray(exterior),
                    'interiors': []
                    }
                e0, e1 = [list(t) for t in zip(*exterior)]
                path = Path(vertices[e0 + [e0[0]], :], closed=True)
        return _index_ring_collection

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
    def logger(self):
        return logging.getLogger(__name__ + '.' + self.__class__.__name__)

    @description.setter
    def description(self, description):
        self._description = description

    def _certify_input_geom(self, geom_type, geom):
        geom_types = {
            "triangles": 3,
            "quads": 4,
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
        elements_id = range(1, len(elements_id)+1) \
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
            "nodes": self._nodes,
            "elements": self._elements,
            "description": description,
        }

    @property
    @lru_cache(maxsize=None)
    def _sms2dm(self):
        description = self.description
        if self.crs is not None:
            description += f"; {self.crs.srs}"
        return {
            "ND": self._nodes,
            "E3T": self._triangles,
            "E4Q": self._quads,
        }

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
        type(self).values.fget.cache_clear()

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
