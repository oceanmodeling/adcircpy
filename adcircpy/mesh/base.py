from abc import ABC
from collections import defaultdict
from functools import lru_cache
from itertools import permutations
import logging
import os
import pathlib
from typing import Dict, Hashable, List, Sequence, Union

import geopandas as gpd
from matplotlib.collections import PolyCollection
from matplotlib.path import Path
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from matplotlib.tri import Triangulation
import numpy as np
from pyproj import CRS, Transformer
from shapely.geometry import box, LinearRing, LineString, MultiPolygon, Point, Polygon

from adcircpy.figures import figure
from adcircpy.mesh.parsers import grd, sms2dm

_logger = logging.getLogger(__name__)


class Nodes:
    def __init__(self, nodes: Dict[Hashable, List[List]], crs=None):
        """Setter for the nodes attribute.

        Argument nodes must be of the form:
            {id: [(x0, y0), z0]}
            or
            {id: [(x0, y0), [z0, ..., zn]}

        Grd format is assumed to be exclusively a 2D format that can hold
        triangles or quads.

        """

        for coords, _ in nodes.values():
            if len(coords) != 2:
                raise ValueError(
                    'Coordinate vertices for a gr3 type must be 2D, but got '
                    f'coordinates {coords}.'
                )

        self._id = list(nodes.keys())
        self._coords = np.array([coords for coords, _ in nodes.values()])
        self._crs = CRS.from_user_input(crs) if crs is not None else crs
        self._values = np.array([value for _, value in nodes.values()])

    def transform_to(self, dst_crs):
        dst_crs = CRS.from_user_input(dst_crs)
        if not self.crs.equals(dst_crs):
            self._coords = self.get_xy(dst_crs)
            self._crs = dst_crs

        if hasattr(self, '_gdf'):
            del self._gdf

    def get_xy(self, crs: Union[CRS, str] = None):
        if crs is not None:
            crs = CRS.from_user_input(crs)
            if not crs.equals(self.crs):
                transformer = Transformer.from_crs(self.crs, crs, always_xy=True)
                x, y = transformer.transform(self.coord[:, 0], self.coord[:, 1])
                return np.vstack([x, y]).T
        return self.coord

    @property
    def gdf(self):
        if not hasattr(self, '_gdf'):
            data = []
            for id, coord, values in zip(self._id, self._coords, self.values):
                data.append({'geometry': Point(coord), 'id': id, 'values': values})
            self._gdf = gpd.GeoDataFrame(data, crs=self.crs)
        return self._gdf

    @property
    def id(self):
        return self._id

    @property
    def index(self):
        if not hasattr(self, '_index'):
            self._index = np.arange(len(self._id))
        return self._index

    @property
    def crs(self):
        return self._crs

    @property
    def values(self):
        return self._values

    @property
    def coords(self):
        return self._coords

    @property
    def coord(self):
        return self.coords

    def get_index_by_id(self, id: Hashable):
        if not hasattr(self, 'node_id_to_index'):
            self.node_id_to_index = {self.id[i]: i for i in range(len(self.id))}
        return self.node_id_to_index[id]

    def get_id_by_index(self, index: int):
        if not hasattr(self, 'node_index_to_id'):
            self.node_index_to_id = {i: self.id[i] for i in range(len(self.id))}
        return self.node_index_to_id[index]

    def to_dict(self):
        nodes = {nid: (coo, val) for nid, coo, val in zip(self._id, self._coords, self.values)}
        return nodes


class Elements:
    def __init__(self, nodes: Nodes, elements: Dict[Hashable, Sequence]):
        if not isinstance(elements, dict):
            raise TypeError('Argument elements must be a dict.')

        vertex_id_set = set(nodes.id)
        for id, geom in elements.items():
            if not isinstance(geom, Sequence):
                raise TypeError(
                    f'Element with id {id} of the elements '
                    f'argument must be of type {Sequence}, not '
                    f'type {type(geom)}.'
                )
            if not set(geom).issubset(vertex_id_set):
                ValueError(f'Element with id {id} is not a subset of the ' "coordinate id's.")
        self.nodes = nodes
        self.elements = elements

    @property
    def id(self):
        if not hasattr(self, '_id'):
            self._id = list(self.elements.keys())
        return self._id

    @property
    def index(self):
        if not hasattr(self, '_index'):
            self._index = np.arange(len(self.elements))
        return self._index

    def get_index_by_id(self, id: Hashable):
        if not hasattr(self, 'element_id_to_index'):
            self.element_id_to_index = {self.id[i]: i for i in range(len(self.id))}
        return self.element_id_to_index[id]

    def get_id_by_index(self, index: int):
        if not hasattr(self, 'element_index_to_id'):
            self.element_index_to_id = {i: self.id[i] for i in range(len(self.id))}
        return self.element_index_to_id[index]

    def get_indexes_around_index(self, index):
        if not hasattr(self, 'indexes_around_index'):

            def append_geom(geom):
                for simplex in geom:
                    for i, j in permutations(simplex, 2):
                        indexes_around_index[i].add(j)

            indexes_around_index = defaultdict(set)
            append_geom(self.triangles)
            append_geom(self.quads)
            self.indexes_around_index = indexes_around_index
        return list(self.indexes_around_index[index])

    def get_ball(self, order: int, id=None, index=None):

        if not isinstance(order, int):
            raise TypeError("Argument 'order' must be of type int.")

        if not order >= 0:
            raise TypeError("Argument 'order' must be of greater " 'than zero.')

        if id is None and index is None:
            raise ValueError('Must specify one keyword argument of index or id.')

        if id is not None and index is not None:
            raise ValueError('Must specify only one keyword argument of ' 'index or id.')

        if id is not None:
            index = self.get_index_by_id(id)

        new_neighbors: Union[list, set]
        eidxs = set([index])
        for i in range(order):
            elements = self.array[list(sorted(eidxs)), :]
            new_neighbors = list(
                map(self.get_indexes_around_index, list(set(elements.data.flatten())))
            )
            new_neighbors = set([item for sublist in new_neighbors for item in sublist])
            eidxs = eidxs.union(
                set(
                    np.where(
                        np.logical_and(
                            np.any(np.isin(self.array, list(set(new_neighbors))), axis=1),
                            np.any(np.isin(self.array, elements), axis=1),
                        )
                    )[0]
                )
            )
        return self.gdf.loc[eidxs].geometry.unary_union.exterior

    @property
    def array(self):
        if not hasattr(self, '_array'):
            rank = int(max(map(len, self.elements.values())))
            array = np.full((len(self.elements), rank), -1)
            for i, element in enumerate(self.elements.values()):
                row = np.array(list(map(self.nodes.get_index_by_id, element)))
                array[i, : len(row)] = row
            array = np.ma.masked_equal(array, -1)
            self._array = array
        return self._array

    @property
    def triangles(self):
        if not hasattr(self, '_triangles'):
            self._triangles = np.array(
                [
                    list(map(self.nodes.get_index_by_id, element))
                    for element in self.elements.values()
                    if len(element) == 3
                ]
            )
        return self._triangles

    @property
    def quadrilaterals(self):
        return self.quads

    @property
    def quads(self):
        if not hasattr(self, '_quads'):
            self._quads = np.array(
                [
                    list(map(self.nodes.get_index_by_id, element))
                    for element in self.elements.values()
                    if len(element) == 4
                ]
            )
        return self._quads

    @property
    def triangulation(self):
        if not hasattr(self, '_triangulation'):
            triangles = self.triangles.tolist()
            for quad in self.quads:
                triangles.append([quad[0], quad[1], quad[3]])
                triangles.append([quad[1], quad[2], quad[3]])
            self._triangulation = Triangulation(
                self.nodes.coord[:, 0], self.nodes.coord[:, 1], triangles
            )
        return self._triangulation

    @property
    def gdf(self):
        if not hasattr(self, '_gdf'):
            data = []
            for id, element in self.elements.items():
                data.append(
                    {
                        'geometry': Polygon(
                            self.nodes.coord[list(map(self.get_index_by_id, element))]
                        ),
                        'id': id,
                    }
                )
            self._gdf = gpd.GeoDataFrame(data, crs=self.nodes.crs)
        return self._gdf


class Edges:
    def __init__(self, grd: 'Grd'):
        self._grd = grd

    @lru_cache(maxsize=1)
    def __call__(self) -> gpd.GeoDataFrame:
        data = []
        for ring in self._grd.hull.rings().itertuples():
            coords = ring.geometry.coords
            for i in range(1, len(coords)):
                data.append(
                    {
                        'geometry': LineString([coords[i - 1], coords[i]]),
                        'bnd_id': ring.bnd_id,
                        'type': ring.type,
                    }
                )
        return gpd.GeoDataFrame(data, crs=self._grd.crs)

    def exterior(self):
        return self().loc[self()['type'] == 'exterior']

    def interior(self):
        return self().loc[self()['type'] == 'interior']


class Rings:
    def __init__(self, grd: 'Grd'):
        self._grd = grd

    @lru_cache(maxsize=1)
    def __call__(self) -> gpd.GeoDataFrame:
        tri = self._grd.elements.triangulation
        idxs = np.vstack(list(np.where(tri.neighbors == -1))).T
        boundary_edges = []
        for i, j in idxs:
            boundary_edges.append((tri.triangles[i, j], tri.triangles[i, (j + 1) % 3]))
        sorted_rings = sort_rings(edges_to_rings(boundary_edges), self._grd.nodes.coord)
        data = []
        for bnd_id, rings in sorted_rings.items():
            coords = self._grd.nodes.coord[rings['exterior'][:, 0], :]
            geometry = LinearRing(coords)
            data.append({'geometry': geometry, 'bnd_id': bnd_id, 'type': 'exterior'})
            for interior in rings['interiors']:
                coords = self._grd.nodes.coord[interior[:, 0], :]
                geometry = LinearRing(coords)
                data.append({'geometry': geometry, 'bnd_id': bnd_id, 'type': 'interior'})
        return gpd.GeoDataFrame(data, crs=self._grd.crs)

    def exterior(self):
        return self().loc[self()['type'] == 'exterior']

    def interior(self):
        return self().loc[self()['type'] == 'interior']


class Hull:
    def __init__(self, grd: 'Grd'):
        self._grd = grd
        self.edges = Edges(grd)
        self.rings = Rings(grd)

    @lru_cache(maxsize=1)
    def __call__(self) -> gpd.GeoDataFrame:
        data = []
        for bnd_id in np.unique(self.rings()['bnd_id'].tolist()):
            exterior = self.rings().loc[
                (self.rings()['bnd_id'] == bnd_id) & (self.rings()['type'] == 'exterior')
            ]
            interiors = self.rings().loc[
                (self.rings()['bnd_id'] == bnd_id) & (self.rings()['type'] == 'interior')
            ]
            data.append(
                {
                    'geometry': Polygon(
                        exterior.iloc[0].geometry.coords,
                        [row.geometry.coords for _, row in interiors.iterrows()],
                    ),
                    'bnd_id': bnd_id,
                }
            )
        return gpd.GeoDataFrame(data, crs=self._grd.crs)

    @lru_cache(maxsize=1)
    def exterior(self):
        data = []
        for exterior in self.rings().loc[self.rings()['type'] == 'exterior'].itertuples():
            data.append({'geometry': Polygon(exterior.geometry.coords)})
        return gpd.GeoDataFrame(data, crs=self._grd.crs)

    @lru_cache(maxsize=1)
    def interior(self):
        data = []
        for interior in self.rings().loc[self.rings()['type'] == 'interior'].itertuples():
            data.append({'geometry': Polygon(interior.geometry.coords)})
        return gpd.GeoDataFrame(data, crs=self._grd.crs)

    @lru_cache(maxsize=1)
    def implode(self) -> gpd.GeoDataFrame:
        return gpd.GeoDataFrame(
            {'geometry': MultiPolygon([polygon.geometry for polygon in self().itertuples()])},
            crs=self._grd.crs,
        )

    @lru_cache(maxsize=1)
    def multipolygon(self) -> MultiPolygon:
        geometry = self.implode().iloc[0].geometry
        if not isinstance(geometry, MultiPolygon):
            geometry = MultiPolygon([geometry])
        return geometry


class Grd(ABC):
    def __init__(self, nodes, elements=None, description=None, crs=None):

        self.nodes = Nodes(nodes, crs)
        self.elements = Elements(self.nodes, elements)
        self.description = "" if description is None else str(description)
        self.hull = Hull(self)

    def __str__(self):
        return grd.to_string(**self.to_dict())

    def to_dict(self):
        return {
            'description': self.description,
            'nodes': self.nodes.to_dict(),
            'elements': self.elements.elements,
            'crs': self.crs,
        }

    def write(self, path, overwrite=False, format='gr3'):
        if format in ['gr3', 'grd']:
            grd.write(self.to_dict(), path, overwrite)
        elif format in ['sms', '2dm', 'sms2dm']:
            sms2dm.write(
                {
                    'ND': {
                        i
                        + 1: (
                            coord,
                            -self.values[i] if not np.isnan(self.values[i]) else -99999,
                        )
                        for i, coord in enumerate(self.coords)
                    },
                    'E3T': {i + 1: index + 1 for i, index in enumerate(self.triangles)},
                    'E4Q': {i + 1: index + 1 for i, index in enumerate(self.quads)},
                },
                path,
                overwrite,
            )
        else:
            raise ValueError(f'Unknown format {format} for hgrid output.')

    def get_xy(self, crs: Union[CRS, str] = None):
        return self.nodes.get_xy(crs)

    def get_bbox(
        self, crs: Union[str, CRS] = None, output_type: str = None
    ) -> Union[Polygon, Bbox]:
        output_type = 'polygon' if output_type is None else output_type
        xmin, xmax = np.min(self.coord[:, 0]), np.max(self.coord[:, 0])
        ymin, ymax = np.min(self.coord[:, 1]), np.max(self.coord[:, 1])
        crs = self.crs if crs is None else crs
        if crs is not None:
            if not self.crs.equals(crs):
                transformer = Transformer.from_crs(self.crs, crs, always_xy=True)
                (xmin, xmax), (ymin, ymax) = transformer.transform((xmin, xmax), (ymin, ymax))
        if output_type == 'polygon':
            return box(xmin, ymin, xmax, ymax)
        elif output_type == 'bbox':
            return Bbox([[xmin, ymin], [xmax, ymax]])
        else:
            raise TypeError(
                "Argument output_type must a string literal 'polygon' or " "'bbox'"
            )

    def invert_sign(self):
        self.nodes.values[:] = -self.nodes.values

    def transform_to(self, dst_crs):
        """Transforms coordinate system of mesh in-place."""
        self.nodes.transform_to(dst_crs)

    def vertices_around_vertex(self, index):
        return self.nodes.vertices_around_vertex(index)

    def copy(self):
        return self.__class__(**self.to_dict())

    @classmethod
    def open(cls, file: Union[str, os.PathLike], crs: Union[str, CRS] = None):
        return cls(**grd.read(pathlib.Path(file), boundaries=False))

    @figure
    def tricontourf(self, axes=None, show=True, figsize=None, cbar=False, **kwargs):
        if len(self.triangles) > 0:
            ax = axes.tricontourf(self.x, self.y, self.triangles, self.values, **kwargs)
            if cbar is True:
                plt.colorbar(ax)
        return axes

    @figure
    def tripcolor(self, axes=None, show=True, figsize=None, **kwargs):
        if len(self.triangles) > 0:
            axes.tripcolor(self.x, self.y, self.triangles, self.values, **kwargs)
        return axes

    @figure
    def triplot(
        self, axes=None, show=False, figsize=None, linewidth=0.07, color='black', **kwargs,
    ):
        if len(self.triangles) > 0:
            kwargs.update({'linewidth': linewidth})
            kwargs.update({'color': color})
            axes.triplot(self.x, self.y, self.triangles, **kwargs)
        return axes

    @figure
    def quadplot(
        self,
        axes=None,
        show=False,
        figsize=None,
        facecolor='none',
        edgecolor='k',
        linewidth=0.07,
        **kwargs,
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

    @figure
    def quadface(self, axes=None, show=False, figsize=None, **kwargs):
        if len(self.quads) > 0:
            pc = PolyCollection(self.coords[self.quads], **kwargs)
            quad_value = np.mean(self.values[self.quads], axis=1)
            pc.set_array(quad_value)
            axes.add_collection(pc)
        return axes

    @figure
    def wireframe(self, axes=None, show=False, **kwargs):
        axes = self.triplot(axes=axes, **kwargs)
        axes = self.quadplot(axes=axes, **kwargs)
        return axes

    @property
    def coords(self):
        return self.nodes.coord

    @property
    def coord(self):
        return self.nodes.coord

    @property
    def vertices(self):
        return self.nodes.coord

    @property
    def vertex_id(self):
        return self.nodes.id

    @property
    def element_id(self):
        return self.elements.id

    @property
    def values(self):
        return self.nodes.values

    @property
    def crs(self):
        return self.nodes.crs

    @property
    def x(self):
        return self.nodes.coord[:, 0]

    @property
    def y(self):
        return self.nodes.coord[:, 1]

    @property
    def triangles(self):
        return self.elements.triangles

    @property
    def quads(self):
        return self.elements.quads

    @property
    def triangulation(self):
        return self.elements.triangulation

    @property
    def bbox(self):
        return self.get_bbox()


def edges_to_rings(edges):
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
            ordered_edges.append(list(reversed(edges.pop(idx))))
        elif ordered_edges[0][0] in e0:
            idx = e0.index(ordered_edges[0][0])
            ordered_edges.insert(0, list(reversed(edges.pop(idx))))
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


def sort_rings(index_rings, vertices):
    """Sorts a list of index-rings.

    Takes a list of unsorted index rings and sorts them into an "exterior" and
    "interior" components. Any doubly-nested rings are considered exterior
    rings.

    TODO: Refactor and optimize. Calls that use :class:matplotlib.path.Path can
    probably be optimized using shapely.
    """

    # sort index_rings into corresponding "polygons"
    areas = list()
    for index_ring in index_rings:
        e0, e1 = [list(t) for t in zip(*index_ring)]
        areas.append(float(Polygon(vertices[e0, :]).area))

    # maximum area must be main mesh
    idx = areas.index(np.max(areas))
    exterior = index_rings.pop(idx)
    areas.pop(idx)
    _id = 0
    _index_rings = dict()
    _index_rings[_id] = {'exterior': np.asarray(exterior), 'interiors': []}
    e0, e1 = [list(t) for t in zip(*exterior)]
    path = Path(vertices[e0 + [e0[0]], :], closed=True)
    while len(index_rings) > 0:
        # find all internal rings
        potential_interiors = list()
        for i, index_ring in enumerate(index_rings):
            e0, e1 = [list(t) for t in zip(*index_ring)]
            if path.contains_point(vertices[e0[0], :]):
                potential_interiors.append(i)
        # filter out nested rings
        real_interiors = list()
        for i, p_interior in reversed(list(enumerate(potential_interiors))):
            _p_interior = index_rings[p_interior]
            check = [
                index_rings[k]
                for j, k in reversed(list(enumerate(potential_interiors)))
                if i != j
            ]
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
            _index_rings[_id]['interiors'].append(np.asarray(index_rings.pop(i)))
            areas.pop(i)
        # if no internal rings found, initialize next polygon
        if len(index_rings) > 0:
            idx = areas.index(np.max(areas))
            exterior = index_rings.pop(idx)
            areas.pop(idx)
            _id += 1
            _index_rings[_id] = {'exterior': np.asarray(exterior), 'interiors': []}
            e0, e1 = [list(t) for t in zip(*exterior)]
            path = Path(vertices[e0 + [e0[0]], :], closed=True)
    return _index_rings


def signed_polygon_area(vertices):
    # https://code.activestate.com/recipes/578047-area-of-polygon-using-shoelace-formula/
    n = len(vertices)  # of vertices
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += vertices[i][0] * vertices[j][1]
        area -= vertices[j][0] * vertices[i][1]
        return area / 2.0
