from abc import ABC
from collections import defaultdict
from functools import lru_cache
from itertools import permutations
import logging
import os
import pathlib
from typing import Hashable, Mapping, Union

import geopandas as gpd
from geopandas import GeoDataFrame
from matplotlib.collections import PolyCollection
from matplotlib.path import Path
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from matplotlib.tri import Triangulation
import numpy
import numpy as np
import pandas
from pandas import DataFrame
from pyproj import CRS, Transformer
from shapely.geometry import box, LinearRing, LineString, MultiPolygon, Polygon
from shapely.ops import polygonize

from adcircpy.figures import figure
from adcircpy.mesh.parsers import grd, sms2dm


class Elements:
    def __init__(
        self,
        nodes: DataFrame,
        elements: DataFrame,
        crs: CRS = None,
        check_elements: bool = False,
    ):
        if check_elements:
            logging.debug('validating elements...')
            vertex_id_set = nodes.index
            for id, geom in elements.iterrows():
                if not np.all(np.in1d(geom[1:][~pandas.isna(geom[1:])].values, vertex_id_set)):
                    ValueError(
                        f'element "{id}" references nodes that do not exist ("{geom[1:].values}")'
                    )

        self.nodes = nodes
        self.elements = elements
        self.crs = crs

    @property
    def id(self):
        if not hasattr(self, '_id'):
            self._id = self.elements.index.values
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
            self._array = self.nodes.values
            # rank = int(max(map(len, self.elements.values())))
            # array = np.full((len(self.elements), rank), -1)
            # for i, element in enumerate(self.elements.values()):
            #     row = np.array(list(map(self.nodes.get_index_by_id, element)))
            #     array[i, : len(row)] = row
            # array = np.ma.masked_equal(array, -1)
            # self._array = array
        return self._array

    @property
    def triangles(self):
        if not hasattr(self, '_triangles'):
            num_nodes = 3
            self._triangles = self.elements.loc[self.elements.iloc[:, 0] == num_nodes].iloc[
                :, 1 : num_nodes + 1
            ]
        return self._triangles

    @property
    def quadrilaterals(self):
        return self.quads

    @property
    def quads(self):
        if not hasattr(self, '_quads'):
            num_nodes = 4
            self._quads = self.elements.loc[self.elements.iloc[:, 0] == num_nodes].iloc[
                :, 1 : num_nodes + 1
            ]
        return self._quads

    @property
    def triangulation(self):
        if not hasattr(self, '_triangulation'):
            triangles = self.triangles.values
            quads = []
            for quad in self.quads.values:
                quads.append([quad[0], quad[1], quad[3]])
                quads.append([quad[1], quad[2], quad[3]])
            if len(quads) > 0:
                triangles = numpy.concatenate([triangles, numpy.array(quads)])
            triangles -= 1
            self._triangulation = Triangulation(
                x=self.nodes.iloc[:, 0].values,
                y=self.nodes.iloc[:, 1].values,
                triangles=triangles,
            )
        return self._triangulation

    @property
    def gdf(self):
        if not hasattr(self, '_gdf'):
            data = [
                {
                    'geometry': Polygon(
                        self.nodes.coord.iloc[
                            element[1:][~pandas.isna(element[1:])].astype(int)
                        ]
                    ),
                    'id': id,
                }
                for id, element in self.elements.iterrows()
            ]
            self._gdf = gpd.GeoDataFrame(data, crs=self.crs)
        return self._gdf


class Hull:
    def __init__(self, grd: 'Grd'):
        self._grd = grd
        self.rings = Rings(self._grd)

    @property
    def edges(self) -> GeoDataFrame:
        data = []
        for ring in self.rings.itertuples():
            coords = ring.geometry.coords
            for i in range(1, len(coords)):
                data.append(
                    {
                        'geometry': LineString([coords[i - 1], coords[i]]),
                        'bnd_id': ring.bnd_id,
                        'type': ring.type,
                    }
                )
        return GeoDataFrame(data, crs=self._grd.crs)


class Rings:
    def __init__(self, grd: 'Grd'):
        self._grd = grd

    @lru_cache(maxsize=1)
    def __call__(self) -> gpd.GeoDataFrame:
        data = []
        for bnd_id, rings in self.sorted().items():
            geometry = LinearRing(self._grd.coords.iloc[rings['exterior'][:, 0], :].values)
            data.append({'geometry': geometry, 'bnd_id': bnd_id, 'type': 'exterior'})
            for interior in rings['interiors']:
                geometry = LinearRing(self._grd.coords.iloc[interior[:, 0], :].values)
                data.append({'geometry': geometry, 'bnd_id': bnd_id, 'type': 'interior'})
        return GeoDataFrame(data, crs=self._grd.crs)

    def exterior(self):
        return self().loc[self()['type'] == 'exterior']

    def interior(self):
        return self().loc[self()['type'] == 'interior']

    @lru_cache(maxsize=1)
    def sorted(self):
        tri = self._grd.elements.triangulation
        idxs = np.vstack(list(np.where(tri.neighbors == -1))).T
        boundary_edges = []
        for i, j in idxs:
            boundary_edges.append((tri.triangles[i, j], tri.triangles[i, (j + 1) % 3]))
        return sort_rings(edges_to_rings(boundary_edges), self._grd.coord.values)

    @property
    def geodataframe(self) -> GeoDataFrame:
        data = []
        rings = self.rings
        for bnd_id in np.unique(rings['bnd_id'].tolist()):
            exterior = rings.loc[(rings['bnd_id'] == bnd_id) & (rings['type'] == 'exterior')]
            interiors = rings.loc[(rings['bnd_id'] == bnd_id) & (rings['type'] == 'interior')]
            data.append(
                {
                    'geometry': Polygon(
                        exterior.iloc[0].geometry.coords,
                        [row.geometry.coords for _, row in interiors.iterrows()],
                    ),
                    'bnd_id': bnd_id,
                }
            )
        return GeoDataFrame(data, crs=self._grd.crs)

    @property
    def exterior(self) -> GeoDataFrame:
        data = []
        for exterior in self.rings().loc[self.rings()['type'] == 'exterior'].itertuples():
            data.append({'geometry': Polygon(exterior.geometry.coords)})
        return GeoDataFrame(data, crs=self._grd.crs)

    @property
    def interior(self) -> GeoDataFrame:
        data = []
        for interior in self.rings().loc[self.rings()['type'] == 'interior'].itertuples():
            data.append({'geometry': Polygon(interior.geometry.coords)})
        return GeoDataFrame(data, crs=self._grd.crs)

    @property
    def implode(self) -> GeoDataFrame:
        return GeoDataFrame(
            {
                'geometry': MultiPolygon(
                    [polygon.geometry for polygon in self.geodataframe.itertuples()]
                )
            },
            crs=self._grd.crs,
        )

    @property
    @lru_cache(maxsize=1)
    def multipolygon(self) -> MultiPolygon:
        triangles = self._grd.elements.triangulation.triangles

        triangle_edges, counts = numpy.unique(
            numpy.sort(
                numpy.concatenate(
                    [triangles[:, :2], triangles[:, 1:], triangles[:, [0, 2]]], axis=0
                ),
                axis=1,
            ),
            axis=0,
            return_counts=True,
        )
        boundary_edges = triangle_edges[counts == 1]

        boundary_edge_points = self._grd.nodes.iloc[:, :2].values[boundary_edges]

        exterior_polygons = collect_interiors(list(polygonize(boundary_edge_points.tolist())))

        coords = self._grd.nodes.values
        x = coords[:, 0]
        y = coords[:, 1]
        total_triangle_area = numpy.sum(
            numpy.abs(
                (
                    x[triangles[:, 0]] * (y[triangles[:, 1]] - y[triangles[:, 2]])
                    + x[triangles[:, 1]] * (y[triangles[:, 2]] - y[triangles[:, 0]])
                    + x[triangles[:, 2]] * (y[triangles[:, 0]] - y[triangles[:, 1]])
                )
                / 2
            )
        )

        if not numpy.isclose(exterior_polygons[-1].area, total_triangle_area):
            polygon_collection = []
            coords = self._grd.coords.values
            for rings in self.sorted().values():
                exterior = coords[rings['exterior'][:, 0], :]
                interiors = []
                for interior in rings['interiors']:
                    interiors.append(coords[interior[:, 0], :])
                polygon_collection.append(Polygon(exterior, interiors))

            exterior_polygons.extend(polygon_collection)
            exterior_polygons = collect_interiors(exterior_polygons)

        multipolygon = MultiPolygon(exterior_polygons)
        if not multipolygon.is_valid:
            try:
                multipolygon = multipolygon.buffer(0)
            except Exception as error:
                logging.exception(error)

        return multipolygon


class Grd(ABC):
    def __init__(
        self,
        nodes: DataFrame,
        elements: DataFrame = None,
        description: str = None,
        crs: CRS = None,
    ):
        if isinstance(nodes, Mapping):
            nodes = DataFrame.from_dict(nodes, orient='index')
        if isinstance(elements, Mapping):
            elements = DataFrame.from_dict(elements, orient='index')

        if crs is not None and not isinstance(crs, CRS):
            crs = CRS.from_user_input(crs)

        records_with_extra_values = np.any(~pandas.isna(nodes.iloc[:, 3:]), axis=1)
        if np.any(records_with_extra_values):
            raise ValueError(
                f'Coordinate vertices for a gr3 type must be 2D, but got coordinates {nodes[records_with_extra_values]}.'
            )

        nodes = nodes.loc[
            :,
            nodes.columns[:2].to_list()
            + nodes.columns[2:][np.any(~pandas.isna(nodes.iloc[:, 2:]), axis=0)].to_list(),
        ]

        self._coords = nodes.iloc[:, :2]
        self._crs = crs
        self._values = nodes.iloc[:, 2:]

        self.nodes = nodes
        self.elements = Elements(self.nodes, elements, crs)
        self.description = "" if description is None else str(description)
        self.hull = Hull(self)

    def __str__(self):
        return grd.to_string(**self.to_dict())

    def to_dict(self):
        return {
            'description': self.description,
            'nodes': self.nodes,
            'elements': self.elements.elements,
            'crs': self.crs,
        }

    def write(self, path, overwrite=False, format='gr3'):
        if format in ['gr3', 'grd']:
            grd.write(self.to_dict(), path, overwrite)
        elif format in ['sms', '2dm', 'sms2dm']:
            nodes = self.nodes.copy()
            nodes[pandas.isna(nodes)] = -99999
            sms2dm.write(
                {'ND': nodes, 'E3T': self.triangles, 'E4Q': self.quads}, path, overwrite
            )
        else:
            raise ValueError(f'Unknown format {format} for hgrid output.')

    def get_xy(self, crs: Union[CRS, str] = None):
        projected_coordinates = self.coords
        if crs is not None:
            crs = CRS.from_user_input(crs)
            if not crs.equals(self.crs):
                transformer = Transformer.from_crs(self.crs, crs, always_xy=True)
                x, y = transformer.transform(
                    projected_coordinates.iloc[:, 0].values,
                    projected_coordinates.iloc[:, 1].values,
                )
                projected_coordinates = np.vstack([x, y]).T
        return projected_coordinates

    def get_bbox(
        self, crs: Union[str, CRS] = None, output_type: str = None
    ) -> Union[Polygon, Bbox]:
        output_type = 'polygon' if output_type is None else output_type
        xmin, xmax = np.min(self.coords.iloc[:, 0]), np.max(self.coords.iloc[:, 0])
        ymin, ymax = np.min(self.coords.iloc[:, 1]), np.max(self.coords.iloc[:, 1])
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

    def transform_to(self, dst_crs: CRS):
        """Transforms coordinate system of mesh in-place."""

        dst_crs = CRS.from_user_input(dst_crs)
        if not self.crs.equals(dst_crs):
            self._coords = self.get_xy(dst_crs)
            self._crs = dst_crs

        if hasattr(self, '_gdf'):
            del self._gdf

    def copy(self):
        return self.__copy__()

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
        return self._coords

    @property
    def coord(self):
        return self.coords

    @property
    def vertices(self):
        return self.coords

    @property
    def vertex_id(self):
        return self.nodes.index

    @property
    def element_id(self):
        return self.elements.index

    @property
    def values(self):
        return self._values

    @property
    def crs(self):
        return self._crs

    @property
    def x(self):
        return self.coords.iloc[:, 0]

    @property
    def y(self):
        return self.coords.iloc[:, 1]

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

    def __copy__(self) -> 'Grd':
        return self.__class__(**self.to_dict())

    def __eq__(self, other: 'Grd') -> bool:
        return self.nodes.equals(other.nodes)


def edges_to_rings(edges: [(int, int)]) -> numpy.ndarray:
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


def sort_rings(index_rings: [[int]], vertices: numpy.ndarray):
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
        if len(index_ring) > 2:
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


def collect_interiors(polygons: [Polygon]) -> [Polygon]:
    multipolygon_indices = []
    for polygon_index in range(len(polygons)):
        if isinstance(polygons[polygon_index], MultiPolygon):
            polygons.extend(polygons[polygon_index])
            multipolygon_indices.append(polygon_index)
    for index in reversed(multipolygon_indices):
        polygons.remove(index)

    polygon_areas = {
        polygon_index: polygon.area for polygon_index, polygon in enumerate(polygons)
    }
    polygons = [polygons[index] for index in sorted(polygon_areas, key=polygon_areas.get)]

    containers = {polygon_index: None for polygon_index in range(len(polygons))}
    for polygon_index, polygon in enumerate(polygons):
        container_index = containers[polygon_index]
        for other_index, other in (
            (other_index, other)
            for other_index, other in enumerate(polygons)
            if other_index != polygon_index and other_index != container_index
        ):
            # find the smallest polygon containing this
            if polygon.within(other):
                if container_index is None or other.area < polygons[container_index].area:
                    container_index = other_index
        containers[polygon_index] = container_index

    return create_polygons(
        hierarchy=container_hierarchy(containers=containers), polygons=polygons,
    )


def container_hierarchy(
    containers: {int: int}, global_container_index: int = None,
) -> {int: {}}:
    hierarchy = {}

    for interior_index, container_index in containers.items():
        if container_index is None:
            hierarchy[interior_index] = container_hierarchy(
                containers={
                    k: v if v != interior_index else None
                    for k, v in containers.items()
                    if v is not None
                },
                global_container_index=interior_index,
            )
        elif container_index == global_container_index:
            hierarchy[container_index] = interior_index

    return hierarchy


def create_polygons(hierarchy: {int: {}}, polygons: [Polygon]) -> [Polygon]:
    exteriors = []
    for exterior_index, interior_indices in hierarchy.items():
        # construct polygon from exterior and holes
        exteriors.append(
            Polygon(
                polygons[exterior_index].exterior.coords,
                holes=list(polygons[exterior_index].interiors)
                + [
                    polygons[interior_index].exterior.coords
                    for interior_index in interior_indices
                ],
            )
        )

        # invert interior holes into islands
        for internal_exterior_indices in interior_indices.values():
            exteriors.extend(
                create_polygons(hierarchy=internal_exterior_indices, polygons=polygons)
            )

    return exteriors
