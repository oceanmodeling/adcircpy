from collections import defaultdict
from functools import lru_cache
import logging
import pathlib

import fiona
from matplotlib import pyplot
from matplotlib.cm import ScalarMappable
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from shapely.geometry import LineString, mapping

from adcircpy.forcing.bctypes import BoundaryCondition
from adcircpy.forcing.winds.base import WindForcing
from adcircpy.mesh import figures as fig, grd, sms2dm
from adcircpy.mesh.base import EuclideanMesh2D


class AdcircMesh(EuclideanMesh2D):
    """
    Class that represents the unstructured planar mesh used by ADCIRC.
    """

    def __init__(
          self,
          nodes,
          elements,
          boundaries=None,
          crs=None,
          description=None,
    ):
        self._nodes = nodes
        self._elements = elements
        super().__init__(**grd.euclidean_mesh({
            'nodes'      : self._nodes,
            'elements'   : self._elements,
            'description': description,
            'crs'        : crs
        }))
        self._boundaries = boundaries

    @staticmethod
    def open(path, crs=None):
        msh = grd.reader(path)
        msh.update({"crs": crs})
        return AdcircMesh(**msh)

    def add_boundary_type(self, ibtype):
        if ibtype not in self.boundaries:
            self.__boundaries[ibtype] = defaultdict()
        else:
            raise Exception(f'Cannot add boundary_type={ibtype}: '
                            f'boundary type already exists.')

    def set_boundary_data(self, ibtype, id, indexes, **properties):
        msg = "Indexes must be subset of node id's."
        for idx in np.asarray(indexes).flatten():
            assert idx in self.node_index.keys(), msg
        self.__boundaries[ibtype][id] = {
            'indexes': indexes,
            **properties
        }

    def delete_boundary_type(self, ibtype):
        del self.__boundaries[ibtype]

    def delete_boundary_data(self, ibtype, id):
        del self.__boundaries[ibtype][id]

    def write_boundaries(self, path, overwrite=False):
        path = pathlib.Path(path)
        if path.exists() and not overwrite:
            msg = "Destination path exists and overwrite=False"
            raise IOError(msg)
        with fiona.open(
              path.absolute(),
              'w',
              driver='ESRI Shapefile',
              crs=self.crs.srs,
              schema={
                  'geometry'  : 'LineString',
                  'properties': {
                      'id'    : 'int',
                      'ibtype': 'str',
                      'bnd_id': 'str'
                  }}) as dst:
            _cnt = 0
            for ibtype, bnds in self.boundaries.items():
                for id, bnd in bnds.items():
                    idxs = list(map(self.get_node_index, bnd['indexes']))
                    linear_ring = LineString(self.xy[idxs].tolist())
                    dst.write({
                        "geometry"  : mapping(linear_ring),
                        "properties": {
                            "id"    : _cnt,
                            "ibtype": ibtype,
                            "bnd_id": f"{ibtype}:{id}"
                        }
                    })
                    _cnt += 1

    def generate_boundaries(
          self,
          threshold=0.,
          land_ibtype=0,
          interior_ibtype=1,
    ):
        if np.any(np.isnan(self.values)):
            raise Exception("Mesh contains invalid values. Raster values must "
                            "be interpolated to the mesh before generating "
                            "boundaries.")

        # generate exterior boundaries
        for ring in self.outer_ring_collection.values():
            # find boundary edges
            edge_tag = np.full(ring.shape, 0)
            edge_tag[np.where(self.values[ring[:, 0]] < threshold)[0], 0] = -1
            edge_tag[np.where(self.values[ring[:, 1]] < threshold)[0], 1] = -1
            edge_tag[np.where(self.values[ring[:, 0]] >= threshold)[0], 0] = 1
            edge_tag[np.where(self.values[ring[:, 1]] >= threshold)[0], 1] = 1
            # sort boundary edges
            ocean_boundary = list()
            land_boundary = list()
            for i, (e0, e1) in enumerate(edge_tag):
                if np.any(np.asarray((e0, e1)) == -1):
                    ocean_boundary.append(tuple(ring[i, :]))
                elif np.any(np.asarray((e0, e1)) == 1):
                    land_boundary.append(tuple(ring[i, :]))
            ocean_boundaries = self.sort_edges(ocean_boundary)
            land_boundaries = self.sort_edges(land_boundary)
            _bnd_id = len(self.boundaries[None])
            for bnd in ocean_boundaries:
                e0, e1 = [list(t) for t in zip(*bnd)]
                e0 = list(map(self.get_node_id, e0))
                data = e0 + [self.get_node_id(e1[-1])]
                self.set_boundary_data(None, _bnd_id, data)
                _bnd_id += 1
            # add land boundaries
            if land_ibtype not in self.boundaries:
                self.add_boundary_type(land_ibtype)
            _bnd_id = len(self._boundaries[land_ibtype])
            for bnd in land_boundaries:
                e0, e1 = [list(t) for t in zip(*bnd)]
                e0 = list(map(self.get_node_id, e0))
                data = e0 + [self.get_node_id(e1[-1])]
                self.set_boundary_data(land_ibtype, _bnd_id, data)
                _bnd_id += 1
        # generate interior boundaries
        _bnd_id = 0
        _interior_boundaries = defaultdict()
        for interiors in self.inner_ring_collection.values():
            for interior in interiors:
                e0, e1 = [list(t) for t in zip(*interior)]
                if self.signed_polygon_area(self.coords[e0, :]) < 0:
                    e0 = list(reversed(e0))
                    e1 = list(reversed(e1))
                e0 = list(map(self.get_node_id, e0))
                e0 += [e0[0]]
                _interior_boundaries[_bnd_id] = e0
                _bnd_id += 1
        self.add_boundary_type(interior_ibtype)
        for bnd_id, data in _interior_boundaries.items():
            self.set_boundary_data(interior_ibtype, bnd_id, data)

    def add_forcing(self, forcing, id=None):
        if isinstance(forcing, BoundaryCondition):
            if id is None:
                for i in range(len(self.open_boundaries)):
                    self.add_forcing(forcing, i)
            else:
                self._boundary_forcing[forcing.btype]["bnd_ids"].add(id)
                self._boundary_forcing[forcing.btype].update({"obj": forcing})
        elif isinstance(forcing, WindForcing):
            self._surface_forcing.update({'imetype': forcing})
        else:
            msg = f"Unrecognized forcing type {forcing}."
            raise Exception(msg)

    # def get_forcing(self, btype):
    #     msg = f"btype must be one of " \
    #           f"{' '.join(list(self._boundary_forcing.keys()))}"
    #     assert btype in
    #     return self._boundary_forcing[btype]

    def add_nodal_attribute(self, name: str, units: str):
        if name in self.get_nodal_attribute_names():
            raise AttributeError(f'Cannot add nodal attribute with name '
                                 f'{name}: attribute already exists.')
        else:
            self.add_attribute(
                name, units=units, coldstart=False, hotstart=False)
            self._nodal_attribute_collection[name] = None

    def set_nodal_attribute(
          self,
          attribute_name,
          values,
          coldstart: bool = False,
          hotstart: bool = False
    ):
        if attribute_name not in self.get_nodal_attribute_names():
            raise AttributeError(f'Cannot set nodal attribute with name '
                                 f'{attribute_name}: attribute has not been '
                                 f'added yet.')
        assert isinstance(coldstart, bool)
        assert isinstance(hotstart, bool)
        properties = {
            'units'    : self._attributes[attribute_name]['units'],
            'coldstart': coldstart,
            'hotstart' : hotstart
        }
        self.set_attribute(attribute_name, values, **properties)

    def get_coldstart_attributes(self):
        coldstart_attributes = dict()
        for attribute in self.get_nodal_attribute_names():
            attr = self.get_attribute(attribute)
            if attr['coldstart']:
                coldstart_attributes[attribute] = attr
        return coldstart_attributes

    def get_hotstart_attributes(self):
        hotstart_attributes = dict()
        for attribute in self.get_nodal_attribute_names():
            attr = self.get_attribute(attribute)
            if attr['hotstart']:
                hotstart_attributes[attribute] = attr
        return hotstart_attributes

    def set_nodal_attribute_coldstart_state(self, attribute, state):
        assert isinstance(state, bool)
        self.get_attribute(attribute)['coldstart'] = state

    def set_nodal_attribute_hotstart_state(self, attribute, state):
        assert isinstance(state, bool)
        self.get_attribute(attribute)['hotstart'] = state

    def set_nodal_attribute_state(
          self,
          attribute,
          coldstart,
          hotstart
    ):
        self.set_nodal_attribute_coldstart_state(attribute, coldstart)
        self.set_nodal_attribute_hotstart_state(attribute, hotstart)

    def get_nodal_attribute_names(self):
        return self._nodal_attribute_collection.keys()

    def get_nodal_attribute(self, name):
        if name not in self.get_nodal_attribute_names():
            msg = f"Nodal attrbiute with name {name} has not been loaded."
            raise AttributeError(msg)
        if self._nodal_attribute_collection[name] is None:
            # TODO: the 'values' array can be generated more succintly.
            def mode_rows(a):
                a = np.ascontiguousarray(a)
                void_dt = np.dtype((np.void,
                                    a.dtype.itemsize * np.prod(a.shape[1:])))
                _, ids, count = np.unique(a.view(void_dt).ravel(),
                                          return_index=True,
                                          return_counts=True)
                largest_count_id = ids[count.argmax()]
                most_frequent_row = a[largest_count_id]
                return most_frequent_row

            _attr = self.get_attribute(name).copy()
            if _attr['values'].ndim == 1:  # rewrite as column major
                _attr['values'] = _attr['values'].reshape(
                    (_attr['values'].shape[0], 1))
            _attr['defaults'] = mode_rows(_attr['values'])
            _attr['non_default_indexes'] = np.where(
                (_attr['values'] != _attr['defaults']).all(axis=1))[0]
            self._nodal_attribute_collection[name] = _attr
        return self._nodal_attribute_collection[name]

    def has_nodal_attribute(self, attribute_name, runtype=None):
        """
        if runtype is None, returns True if attribute_name is in any.
        """
        if attribute_name not in self.get_nodal_attribute_names():
            return False
        else:
            assert runtype in [None, 'coldstart', 'hotstart']
            attr = self.get_attribute(attribute_name)
            if runtype in [None, 'coldstart']:
                if attr['coldstart']:
                    return True
            if runtype in [None, 'hotstart']:
                if attr['hotstart']:
                    return True
            return False

    def import_nodal_attributes(self, fort13):
        fort13 = self._parse_fort13(fort13)
        if fort13.pop('NumOfNodes') != len(self.node_id):
            raise Exception('fort.13 file does not match the mesh.')
        self._AGRID = fort13.pop('AGRID')
        for attribute, data in fort13.items():
            values = np.asarray(data['values'])
            if values.ndim == 1:
                values = values.reshape((values.shape[0], 1))
            full_values = np.full(
                (self.values.size,
                 np.asarray(data['default_values']).flatten().size),
                np.nan)
            for i, idx in enumerate(data['indexes']):
                for j, value in enumerate(values[i, :].tolist()):
                    full_values[idx, j] = value
            idxs = np.where(np.isnan(full_values).all(axis=1))[0]
            for idx in idxs:
                for i, value in enumerate(data['default_values']):
                    full_values[idx, i] = value
            # converts from column major to row major, leave it column major.
            # if full_values.shape[1] == 1:
            #     full_values = full_values.flatten()
            self.add_nodal_attribute(attribute, data['units'])
            self.set_nodal_attribute(attribute, full_values)

    def generate_tau0(
          self,
          default_value=0.03,
          threshold_distance=1750.,
          shallow_tau0=0.02,
          deep_tau0=0.005,
          threshold_depth=-10.,
          coldstart=True,
          hotstart=True
    ):
        """
        Reimplementation of tau0_gen.f by Robert Weaver (2008)
        1) computes  distance to each neighboring node
        2) averages all distances to find rep. distance @ each node.
        3) Assigns a tau0 value based on depth and rep. distance.
        Asssumes threshold_distance is given in meters.
        """
        msg = "Cannot compute TAU0 with nan depth values."
        assert not np.any(np.isnan(self.values)), msg
        msg = "Cannot compute TAU0 with no coordinate reference system set."
        assert self.crs is not None, msg
        points = self.get_xy(3395)
        values = np.full(self.values.shape, default_value)
        for k, v in self.node_neighbors.items():
            x0, y0 = points[k]
            distances = list()
            for idx in v:
                x1, y1 = points[idx]
                distances.append(np.sqrt((x0 - x1) ** 2 + (y0 - y1) ** 2))
            distance = np.mean(distances)
            if distance >= threshold_distance:
                if self.values[k] >= threshold_depth:
                    values[k] = shallow_tau0
                else:
                    values[k] = deep_tau0
        if 'primitive_weighting_in_continuity_equation' \
              not in self.get_nodal_attribute_names():
            self.add_nodal_attribute(
                'primitive_weighting_in_continuity_equation',
                'unitless'
            )
        self.set_nodal_attribute(
            'primitive_weighting_in_continuity_equation',
            values
        )

    def write_fort14(self, path, overwrite=False):
        super().write(path, overwrite, 'gr3')

    def write_fort13(self, path, overwrite=False):
        if path is not None:
            path = pathlib.Path(path)
            if path.is_file() and not overwrite:
                msg = 'File exists, pass overwrite=True to allow overwrite.'
                raise Exception(msg)
            else:
                with open(path, 'w') as f:
                    f.write(self.fort13)
        else:
            print(self.fort13)

    def critical_timestep(self, cfl, maxvel=5., g=9.8, method='simple'):
        """
        http://swash.sourceforge.net/online_doc/swashuse/node47.html
        """
        msg = "method keyword must be 'simple' or 'conservative'"
        assert method in ['simple', 'conservative'], msg
        dxdy = len(self.values) * [None]
        for k, v in self.node_distances_meters.items():
            _dxdy = []
            for idx in v:
                _dxdy.append(self.node_distances_meters[k][idx])
            dxdy[k] = np.min(_dxdy)
        if method == 'simple':
            return cfl * np.min(dxdy) / np.abs(maxvel)
        elif method == 'conservative':
            n = cfl * np.asarray(dxdy)
            d = np.sqrt(g * np.abs(self.values)) + np.abs(maxvel)
            return np.min(np.divide(n, d))

    def limgrad(self, dfdx, imax=100, ftol=None, verbose=False, minimize=True):
        if not self.crs.is_geographic:
            original_crs = self.crs
            self.transform_to("EPSG:3395")
        ffun = self.values.copy()
        edges = self.triangulation.edges
        if ftol is None:
            ftol = np.min(ffun) * np.sqrt(np.finfo(float).eps)
        ftol = float(ftol)

        # precompute distances
        distances = np.sqrt(
            (self.x[edges[:, 0]] - self.x[edges[:, 1]]) ** 2
            +
            (self.y[edges[:, 0]] - self.y[edges[:, 1]]) ** 2
        )
        dz = np.abs(ffun[edges[:, 0]] - ffun[edges[:, 1]])

        def get_active_edges_wetdry():
            idxs = np.where(np.sum(np.sign(self.values[edges]), axis=1) == 0)
            _active_eges = np.zeros(edges.shape[0], dtype=bool)
            _active_eges[idxs] = True
            active_edges = edges[np.where(
                np.logical_and(
                    np.divide(dz, distances) > dfdx,
                    _active_eges
                ))]
            return active_edges

        def get_active_edges_traditional():
            return edges[np.divide(dz, distances) > dfdx]

        active_edges = get_active_edges_traditional()
        # active_edges = get_active_edges_wetdry
        cnt = 0
        cnt_table = [len(active_edges)]
        if verbose:
            print(f"iteration: {cnt}, "
                  f"remaining points: {len(active_edges)}")
        _iter = False
        while len(active_edges) > 0:
            cnt += 1
            active_distances = np.sqrt(
                (self.x[active_edges[:, 0]] - self.x[active_edges[:, 1]]) ** 2
                +
                (self.y[active_edges[:, 0]] - self.y[active_edges[:, 1]]) ** 2
            )
            for i, (p0, p1) in enumerate(active_edges):
                z0, z1 = ffun[p0], ffun[p1]
                # push down method
                # if z0 > z1:
                #     ffun[p0] = z0 - dfdx*active_distances[i]
                # elif z0 < z1:
                #     ffun[p1] = z0 - dfdx*active_distances[i]

                # traditional method
                z0, z1 = ffun[p0], ffun[p1]
                if z0 < z1:
                    ffun[p1] = z0 + dfdx * active_distances[i]
                elif z0 > z1:
                    ffun[p1] = z0 - dfdx * active_distances[i]

            dz = np.abs(ffun[edges[:, 0]] - ffun[edges[:, 1]])
            active_edges = get_active_edges_traditional()
            # active_edges = get_active_edges_wetdry()
            cnt_table.append(len(active_edges))
            if verbose:
                print(f"iteration: {cnt}, "
                      f"remaining points: {len(active_edges)}")
            if imax == cnt:
                break
        if cnt == imax:
            if minimize:
                idx = cnt_table.index(min(cnt_table))
                if idx > 0:
                    self.limgrad(
                        dfdx,
                        idx,
                        verbose=verbose,
                        minimize=False
                    )

        if 'original_crs' in locals():
            self.transform_to(original_crs)
        self._values = ffun

    # plotting functions
    @fig._figure
    def make_plot(
          self,
          axes: pyplot.Axes = None,
          vmin: float = None,
          vmax: float = None,
          show: bool = False,
          title: str = None,
          # figsize=rcParams["figure.figsize"],
          extent: (float, float, float, float) = None,
          cbar_label: str = None,
          **kwargs
    ):
        if vmin is None:
            vmin = np.min(self.values)
        if vmax is None:
            vmax = np.max(self.values)
        kwargs.update(**fig.get_topobathy_kwargs(self.values, vmin, vmax))
        kwargs.pop('col_val')
        levels = kwargs.pop('levels')
        if vmin != vmax:
            self.tricontourf(
                axes=axes,
                levels=levels,
                vmin=vmin,
                vmax=vmax,
                **kwargs
            )
        else:
            self.tripcolor(axes=axes, **kwargs)
        self.quadface(axes=axes, **kwargs)
        axes.axis('scaled')
        if extent is not None:
            axes.axis(extent)
        if title is not None:
            axes.set_title(title)
        mappable = ScalarMappable(cmap=kwargs['cmap'])
        mappable.set_array([])
        mappable.set_clim(vmin, vmax)
        divider = make_axes_locatable(axes)
        cax = divider.append_axes("bottom", size="2%", pad=0.5)
        cbar = plt.colorbar(
            mappable,
            cax=cax,
            orientation='horizontal'
        )
        cbar.set_ticks([vmin, vmax])
        cbar.set_ticklabels([np.around(vmin, 2), np.around(vmax, 2)])
        if cbar_label is not None:
            cbar.set_label(cbar_label)
        if show:
            pyplot.show()
        return axes

    @fig._figure
    def plot_boundary(
          self,
          ibtype: str,
          id: int,
          tags: bool = True,
          axes: pyplot.Axes = None,
          show: bool = False,
          figsize: (float, float) = None,
          **kwargs
    ):

        boundary = list(map(
            self.get_node_index,
            self.boundaries[ibtype][id]['indexes']
        ))
        p = axes.plot(self.x[boundary], self.y[boundary], **kwargs)
        if tags:
            axes.text(
                self.x[boundary[len(boundary) // 2]],
                self.y[boundary[len(boundary) // 2]],
                f"ibtype={ibtype}\nid={id}",
                color=p[-1].get_color()
            )
        if show:
            pyplot.show()
        return axes

    @fig._figure
    def plot_boundaries(
          self,
          axes: pyplot.Axes = None,
          show: bool = False,
          figsize: (float, float) = None,
          **kwargs
    ):
        kwargs.update({'axes': axes})
        for ibtype, bnds in self.boundaries.items():
            for id in bnds:
                axes = self.plot_boundary(ibtype, id, **kwargs)
                kwargs.update({'axes': axes})
        if show:
            pyplot.show()
        return kwargs['axes']

    @property
    def fort14(self):
        f = f'{self.description}\n'
        f += f'{self.NE} {self.NP}\n'
        for i in range(self.NP):
            f += f'{i + 1:d} ' \
                 f'{self.x[i]:<.16E} ' \
                 f'{self.y[i]:<.16E} ' \
                 f'{-self.values[i]:<.16E}\n'
        for i in range(self.NE):
            f += f'{i + 1:d} {len(self.elements[i])} ' \
                 f'{self.elements[i, 0] + 1:d} ' \
                 f'{self.elements[i, 1] + 1:d} ' \
                 f'{self.elements[i, 2] + 1:d}\n'
        # ocean boundaries
        f += f'{len(self.ocean_boundaries):d} '
        f += '! total number of ocean boundaries\n'
        # count total number of ocean boundaries
        _sum = np.sum([len(boundary) for boundary in self.ocean_boundaries])
        f += f'{int(_sum):d} ! total number of ocean boundary nodes\n'
        # write ocean boundary indexes
        for i, boundary in enumerate(self.ocean_boundaries):
            f += f'{len(boundary):d} '
            f += '! number of nodes for ocean_boundary_{i}\n'
            for idx in boundary:
                f += f'{idx + 1:d}\n'
        # count remaining boundaries
        num_remaining_boundaries = sum((
            len(self.land_boundaries),
            len(self.inner_boundaries),
            len(self.inflow_boundaries),
            len(self.outflow_boundaries),
            len(self.weir_boundaries),
            len(self.culvert_boundaries)
        ))
        f += f'{num_remaining_boundaries:d} '
        f += '! total number of non-ocean boundaries\n'
        # count total remaining boundary simplices
        num_remaining_boundary_simplices = int(sum((
            np.sum([len(x) for x in self.land_boundaries]),
            np.sum([len(x) for x in self.inner_boundaries]),
            np.sum([len(x) for x in self.inflow_boundaries]),
            np.sum([len(x) for x in self.outflow_boundaries]),
            np.sum([2 * len(x) for x in self.weir_boundaries]),
            np.sum([2 * len(x) for x in self.culvert_boundaries])
        )))

        f += f'{num_remaining_boundary_simplices:d} '
        f += '! Total number of non-ocean boundaries nodes\n'
        # write land boundaries
        for i, boundary in enumerate(self.land_boundaries):
            ibtype = self.get_attribute('land_boundaries')['ibtype'][i]
            f += f'{len(boundary):d} {ibtype:d} '
            f += '! number of nodes and ibtype for land_boundary_{i}\n'
            for idx in boundary:
                f += f'{idx + 1:d}\n'
        # write inner boundaries
        for i, boundary in enumerate(self.inner_boundaries):
            ibtype = self.get_attribute('inner_boundaries')['ibtype'][i]
            f += f'{len(boundary):d} {ibtype:d} '
            f += '! number of nodes and ibtype for inner_boundary_{i}\n'
            for idx in boundary:
                f += f'{idx + 1:d}\n'
        # inflow boundaries
        for i, boundary in enumerate(self.inflow_boundaries):
            ibtype = self.get_attribute('inflow_boundaries')['ibtype'][i]
            f += f'{len(boundary):d} {ibtype:d} '
            f += '! number of nodes and ibtype for inflow_boundary_{i}\n'
            for idx in boundary:
                f += f'{idx + 1:d}\n'
        # outflow boundaries
        for i, boundary in enumerate(self.outflow_boundaries):
            bdata = self.get_attribute('outflow_boundaries')
            ibtype = bdata['ibtype'][i]
            f += f'{len(boundary):d} {ibtype:d} '
            f += '! number of nodes and ibtype for outflow_boundary_{i}\n'
            for j, idx in enumerate(boundary):
                f += f'{idx + 1:d} ' \
                     f'{bdata["barrier_height"][i][j]:G} ' \
                     f'{bdata["supercritical_flow_coefficient"][i][j]:G}\n'
        # weir boundaries
        for i, boundary in enumerate(self.weir_boundaries):
            bdata = self.get_attribute('weir_boundaries')
            ibtype = bdata['ibtype'][i]
            f += f'{len(boundary):d} {ibtype:d} '
            f += '! number of nodes and ibtype for weir_boundary_{i}\n'
            for j, idx in enumerate(boundary):
                f += f'{idx[0] + 1:d} {idx[1] + 1:d} ' \
                     f'{bdata["barrier_height"][i][j]:G} ' \
                     f'{bdata["subcritical_flow_coefficient"][i][j]:G} ' \
                     f'{bdata["supercritical_flow_coefficient"][i][j]:G}\n'
        # culvert boundaries
        for i, boundary in enumerate(self.culvert_boundaries):
            bdata = self.get_attribute('culvert_boundaries')
            ibtype = bdata['ibtype'][i]
            f += f'{len(boundary):d} {ibtype:d} '
            f += '! number of nodes and ibtype for culvert_boundary_{i}\n'
            for j, idx in enumerate(boundary):
                f += f'{idx[0] + 1:d} ' \
                     f'{idx[1] + 1:d} ' \
                     f'{bdata["barrier_height"][i][j]:G} ' \
                     f'{bdata["subcritical_flow_coefficient"][i][j]:G} ' \
                     f'{bdata["supercritical_flow_coefficient"][i][j]:G} ' \
                     f'{bdata["cross_barrier_pipe_height"][i][j]:G} ' \
                     f'{bdata["friction_factor"][i][j]:G} ' \
                     f'{bdata["pipe_diameter"][i][j]:G}\n'
        # f += f"{}\n".format(self.SpatialReference.ExportToWkt())
        return f

    @property
    def fort13(self):
        f = f'{self.description}\n' \
            f'{len(self.node_id)}\n' \
            f'{len(self.get_nodal_attribute_names())}\n'
        for name in self.get_nodal_attribute_names():
            attribute = self.get_nodal_attribute(name)
            f += f'{name}\n' \
                 f'{attribute["units"]}\n' \
                 f'{len(attribute["defaults"])}\n'
            f += ' '.join(f'{n:<.16E}' for n in attribute['defaults']) + '\n'
        for name in self.get_nodal_attribute_names():
            attribute = self.get_nodal_attribute(name)
            f += f'{name}\n' \
                 f'{len(attribute["non_default_indexes"])}\n'
            for idx in attribute['non_default_indexes']:
                f += f'{idx + 1:d} '
                f += ' '.join(f'{val:<.16E}'
                              for val in attribute['values'][idx, :]) + '\n'
        return f

    @property
    def boundaries(self):
        return self._boundaries

    @property
    def primitive_weighting_in_continuity_equation(self):
        return self.get_attribute("primitive_weighting_in_continuity_equation")

    @property
    def surface_submergence_state(self):
        return self.get_attribute("surface_submergence_state")

    @property
    def quadratic_friction_coefficient_at_sea_floor(self):
        return self.get_attribute(
            "quadratic_friction_coefficient_at_sea_floor")

    @property
    def surface_directional_effective_roughness_length(self):
        return self.get_attribute(
            "surface_directional_effective_roughness_length")

    @property
    def surface_canopy_coefficient(self):
        return self.get_attribute("surface_canopy_coefficient")

    @property
    def bridge_pilings_friction_parameters(self):
        return self.get_attribute("bridge_pilings_friction_parameters")

    @property
    def mannings_n_at_sea_floor(self):
        return self.get_attribute("mannings_n_at_sea_floor")

    @mannings_n_at_sea_floor.setter
    def mannings_n_at_sea_floor(self, mannings_n_at_sea_floor: [float]):
        if not self.has_attribute('mannings_n_at_sea_floor'):
            self.add_nodal_attribute('mannings_n_at_sea_floor', 'meters')
        self.set_nodal_attribute('mannings_n_at_sea_floor',
                                 mannings_n_at_sea_floor)

    @property
    def chezy_friction_coefficient_at_sea_floor(self):
        return self.get_attribute("chezy_friction_coefficient_at_sea_floor")

    @property
    def sea_surface_height_above_geoid(self):
        return self.get_attribute("sea_surface_height_above_geoid")

    @property
    def bottom_roughness_length(self):
        return self.get_attribute("bottom_roughness_length")

    @property
    def wave_refraction_in_swan(self):
        return self.get_attribute("wave_refraction_in_swan")

    @property
    def average_horizontal_eddy_viscosity_in_sea_water_wrt_depth(self):
        return self.get_attribute(
            "average_horizontal_eddy_viscosity_in_sea_water_wrt_depth")

    @property
    def elemental_slope_limiter(self):
        return self.get_attribute("elemental_slope_limiter")

    @property
    def advection_state(self):
        return self.get_attribute("advection_state")

    @property
    def initial_river_elevation(self):
        return self.get_attribute("initial_river_elevation")

    @property
    def ics(self):
        if self.crs is None:
            raise Exception('Can\'t determine ics parameter. No projection '
                            'information has been provided for the hgrid.')
        return 2 if self.crs.is_geographic else 1

    @property
    def slam0(self):
        return np.median(self.get_x('EPSG:4326'))

    @property
    def sfea0(self):
        return np.median(self.get_y('EPSG:4326'))

    @property
    def open_boundaries(self):
        open_boundaries = self._open_boundaries.copy()
        for id in self._open_boundaries:

            indexes = list(map(
                self.get_node_index,
                self.boundaries[None][id]['indexes']))

            def iettype(id):
                if id in self._boundary_forcing['iettype']["bnd_ids"]:
                    return self._boundary_forcing['iettype']["obj"].iettype
                return 0

            def ifltype(id):
                if id in self._boundary_forcing['ifltype']["bnd_ids"]:
                    return self._boundary_forcing['ifltype']["obj"].ifltype
                return 0

            def itetype(id):
                if id in self._boundary_forcing['itetype']["bnd_ids"]:
                    return self._boundary_forcing['itetype']["obj"].itetype
                return 0

            def isatype(id):
                if id in self._boundary_forcing['isatype']["bnd_ids"]:
                    return self._boundary_forcing['isatype']["obj"].isatype
                return 0

            def itrtype(id):
                if id in self._boundary_forcing['itrtype']["bnd_ids"]:
                    return self._boundary_forcing['itrtype']["obj"].itrtype
                return 0

            open_boundaries[id].update({
                'indexes': indexes,
                'neta'   : len(indexes),
                'iettype': iettype(id),
                'ifltype': ifltype(id),
                'itetype': itetype(id),
                'isatype': isatype(id),
                'itrtype': itrtype(id)
            })
        return open_boundaries

    @property
    @lru_cache(maxsize=None)
    def _logger(self):
        return logging.getLogger(__name__ + '.' + self.__class__.__name__)

    @property
    def _grd(self):
        grd = super()._grd
        grd.update({'boundaries': self.boundaries})
        return grd

    @property
    def _sms2dm(self):
        _sms2dm = super()._sms2dm
        _sms2dm.update({'nodestrings': sms2dm.nodestrings(self.boundaries)})
        return _sms2dm

    @property
    def _boundaries(self):
        return self.__boundaries

    @_boundaries.setter
    def _boundaries(self, boundaries):
        del self._boundaries
        if boundaries is not None:
            for ibtype, bnds in boundaries.items():
                if ibtype is not None:
                    self.add_boundary_type(ibtype)
                for id, bnd in bnds.items():
                    if 'properties' in bnd.keys():
                        properties = bnd['properties']
                    else:
                        properties = {}
                    self.set_boundary_data(
                        ibtype,
                        id,
                        bnd['indexes'],
                        **properties
                    )

    @_boundaries.deleter
    def _boundaries(self):
        self.__boundaries = {}
        self.__boundaries[None] = {}

    @property
    def _nodes(self):
        return self.__nodes

    @_nodes.setter
    def _nodes(self, nodes):
        assert isinstance(nodes, dict)
        assert all(
            len(node) == 2 and all(type(value) is float for value in node[0])
            for node in nodes.values())
        self.__nodes = nodes

    @property
    def _elements(self):
        return self.__elements

    @_elements.setter
    def _elements(self, elements):
        assert isinstance(elements, dict)
        for element in elements.values():
            assert len(element) >= 3
            for value in element:
                assert value in self._nodes
        self.__elements = elements

    @property
    @lru_cache(maxsize=None)
    def _nodal_attribute_collection(self):
        return {}

    @property
    @lru_cache(maxsize=None)
    def _open_boundaries(self):
        from copy import deepcopy
        return deepcopy(self.__boundaries[None])

    @property
    @lru_cache(maxsize=None)
    def _boundary_forcing(self):
        return {
            'iettype': {"obj": None, "bnd_ids": set()},  # eta
            'ifltype': {"obj": None, "bnd_ids": set()},  # flow
            'itetype': {"obj": None, "bnd_ids": set()},  # temp
            'isatype': {"obj": None, "bnd_ids": set()},  # salinity
            'itrtype': {"obj": None, "bnd_ids": set()},  # tracer
            'iwrtype': {"obj": None, "bnd_ids": set()},  # waves
        }

    @property
    @lru_cache(maxsize=None)
    def _surface_forcing(self):
        return {
            'imetype': None,  # meteorological,
            'iwrtype': None,
            'isrtype': None,  # rain
        }
