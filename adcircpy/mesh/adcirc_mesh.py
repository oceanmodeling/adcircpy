import pathlib
import uuid
import numpy as np
import matplotlib.pyplot as plt
from haversine import haversine, Unit
# import warnings
# from matplotlib.path import Path
from functools import lru_cache
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
from adcircpy.figures import FixPointNormalize
from adcircpy.mesh.unstructured_mesh import UnstructuredMesh


class AdcircMesh(UnstructuredMesh):

    def __init__(
        self,
        vertices,
        elements,
        values,
        crs=None,
        description=None,
        node_id=None,
        element_id=None,
        ocean_boundaries=None,
        land_boundaries=None,
        inner_boundaries=None,
        inflow_boundaries=None,
        outflow_boundaries=None,
        weir_boundaries=None,
        culvert_boundaries=None,
        fort13=None
    ):
        super().__init__(vertices, elements, crs)
        self._description = description
        self._node_id = node_id
        self._element_id = element_id
        self._values = values
        self._ocean_boundaries = ocean_boundaries
        self._land_boundaries = land_boundaries
        self._inner_boundaries = inner_boundaries
        self._inflow_boundaries = inflow_boundaries
        self._outflow_boundaries = outflow_boundaries
        self._weir_boundaries = weir_boundaries
        self._culvert_boundaries = culvert_boundaries
        self._fort13 = fort13
        # self.pad_lakes(H0=0.05, factor=2.)

    @classmethod
    def open(cls, file, crs=None, fort13=None):
        return cls(**cls.parse_fort14(file), crs=crs, fort13=fort13)

    def add_nodal_attribute(self, name, units):
        if name in self.get_nodal_attribute_names():
            raise AttributeError(
                f'Cannot add nodal attribute with name {name}:'
                + ' attribute already exists.')
        else:
            self.add_attribute(
                name, units=units, coldstart=False, hotstart=False)
            self._nodal_attribute_collection[name] = None

    def set_nodal_attribute(
        self,
        attribute_name,
        values,
        coldstart=False,
        hotstart=False
    ):
        if attribute_name not in self.get_nodal_attribute_names():
            raise AttributeError(
                f'Cannot set nodal attribute with name {attribute_name}:'
                + ' attribute has not been added yet.')
        assert isinstance(coldstart, bool)
        assert isinstance(hotstart, bool)
        properties = {
            'units': self._attributes[attribute_name]['units'],
            'coldstart': coldstart,
            'hotstart': hotstart
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
        return self.nodal_attribute_collection.keys()

    def get_nodal_attribute(self, name):
        if name not in self.get_nodal_attribute_names():
            msg = f"Nodal attrbiute with name {name} has not been loaded."
            raise AttributeError(msg)
        if self.nodal_attribute_collection[name] is None:
            # TODO: the 'values' array can be generated more succintly.
            def mode_rows(a):
                a = np.ascontiguousarray(a)
                void_dt = np.dtype(
                    (np.void, a.dtype.itemsize * np.prod(a.shape[1:])))
                _, ids, count = np.unique(
                    a.view(void_dt).ravel(), return_index=1, return_counts=1)
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
        return self.nodal_attribute_collection[name]

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
        fort13 = self.parse_fort13(fort13)
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
        hotstart=True,
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
                distances.append(np.sqrt((x0 - x1)**2 + (y0 - y1)**2))
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

    def critical_timestep(self, cfl, maxvel=5., g=9.8):
        """
        http://swash.sourceforge.net/online_doc/swashuse/node47.html
        """
        dxdy = len(self.values)*[None]
        for k, v in self.node_distances_meters.items():
            _dxdy = []
            for idx in v:
                _dxdy.append(self.node_distances_meters[k][idx])
            dxdy[k] = np.min(_dxdy)
        n = cfl * np.asarray(dxdy)
        d = np.sqrt(g*np.abs(self.values)) + np.abs(maxvel)
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
            (self.x[edges[:, 0]] - self.x[edges[:, 1]])**2
            +
            (self.y[edges[:, 0]] - self.y[edges[:, 1]])**2
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
            msg = f"iteration: {cnt}, "
            msg += f"remaining points: {len(active_edges)}"
            print(msg)
        _iter = False
        while len(active_edges) > 0:
            cnt += 1
            active_distances = np.sqrt(
                (self.x[active_edges[:, 0]] - self.x[active_edges[:, 1]])**2
                +
                (self.y[active_edges[:, 0]] - self.y[active_edges[:, 1]])**2
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
                    ffun[p1] = z0 + dfdx*active_distances[i]
                elif z0 > z1:
                    ffun[p1] = z0 - dfdx*active_distances[i]

            dz = np.abs(ffun[edges[:, 0]] - ffun[edges[:, 1]])
            active_edges = get_active_edges_traditional()
            # active_edges = get_active_edges_wetdry()
            cnt_table.append(len(active_edges))
            if verbose:
                msg = f"iteration: {cnt}, "
                msg += f"remaining points: {len(active_edges)}"
                print(msg)
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

    def triplot(
        self,
        axes=None,
        show=False,
        linewidth=0.07,
        color='black',
        alpha=0.5,
        **kwargs
    ):
        if axes is None:
            fig = plt.figure()
            axes = fig.add_subplot(111)
        axes.triplot(
            self.triangulation,
            linewidth=linewidth,
            color=color,
            alpha=alpha,
            **kwargs
            )
        if show:
            axes.axis('scaled')
            plt.show()
        return axes

    def make_plot(
        self,
        axes=None,
        vmin=None,
        vmax=None,
        cmap='topobathy',
        levels=None,
        show=False,
        title=None,
        figsize=None,
        colors=256,
        extent=None,
        cbar_label=None,
        norm=None,
        **kwargs
    ):
        if axes is None:
            fig = plt.figure(figsize=figsize)
            axes = fig.add_subplot(111)
        if vmin is None:
            vmin = np.min(self.values)
        if vmax is None:
            vmax = np.max(self.values)
        cmap, norm, levels, col_val = self._get_cmap(
            vmin, vmax, cmap, levels, colors, norm)
        axes.tricontourf(
            self.triangulation,
            self.values,
            levels=levels,
            cmap=cmap,
            norm=norm,
            vmin=vmin,
            vmax=vmax,
            **kwargs
            )
        axes.axis('scaled')
        if extent is not None:
            axes.axis(extent)
        if title is not None:
            axes.set_title(title)
        mappable = ScalarMappable(cmap=cmap)
        mappable.set_array([])
        mappable.set_clim(vmin, vmax)
        divider = make_axes_locatable(axes)
        cax = divider.append_axes("bottom", size="2%", pad=0.5)
        cbar = plt.colorbar(
            mappable,
            cax=cax,
            # extend=cmap_extend,
            orientation='horizontal'
        )
        if col_val != 0:
            cbar.set_ticks([vmin, vmin + col_val * (vmax-vmin), vmax])
            cbar.set_ticklabels([np.around(vmin, 2), 0.0, np.around(vmax, 2)])
        else:
            cbar.set_ticks([vmin, vmax])
            cbar.set_ticklabels([np.around(vmin, 2), np.around(vmax, 2)])
        if cbar_label is not None:
            cbar.set_label(cbar_label)
        # plt.gcf().canvas.mpl_connect(
        #     'button_press_event',
        #     lambda e: print(self.get_value(e.xdata, e.ydata))
        #     )
        if show is True:
            plt.show()
        return axes

    def plot_ocean_boundaries(
        self,
        axes=None,
        show=False,
        title=None,
        **kwargs
    ):
        kwargs.update({
            "axes": axes,
            "show": show,
            "title": title
            })
        self._make_boundary_plot(self.ocean_boundaries, **kwargs)

    def plot_land_boundary(
        self,
        axes=None,
        show=False,
        title=None,
        **kwargs
    ):
        kwargs.update({
            "axes": axes,
            "show": show,
            "title": title
            })
        self._make_boundary_plot(self.land_boundaries, **kwargs)

    def plot_nodal_attribute(self, attribute_name):
        raise NotImplementedError

    def write_fort14(self, path, overwrite=False):
        if path is not None:
            path = pathlib.Path(path)
            if path.is_file() and not overwrite:
                raise Exception(
                    'File exists, pass overwrite=True to allow overwrite.')
            else:
                with open(path, 'w') as f:
                    f.write(self.fort14)
        else:
            print(self.fort14)

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

    def _make_boundary_plot(
        self,
        boundary_collection,
        axes=None,
        show=False,
        title=None,
        **kwargs
    ):
        if axes is None:
            fig = plt.figure()
            axes = fig.add_subplot(111)
        for boundary in boundary_collection:
            axes.plot(self.x[boundary], self.y[boundary], **kwargs)
        if show:
            plt.show()

    def _get_cmap(
        self,
        vmin,
        vmax,
        cmap=None,
        levels=None,
        colors=256,
        norm=None
    ):
        colors = int(colors)
        if cmap is None:
            cmap = plt.cm.get_cmap('jet')
            if levels is None:
                levels = np.linspace(vmin, vmax, colors)
            col_val = 0.
        elif cmap == 'topobathy':
            if vmax <= 0.:
                cmap = plt.cm.seismic
                col_val = 0.
                levels = np.linspace(vmin, vmax, colors)
            else:
                wet_count = int(np.floor(colors*(float((self.values < 0.).sum())
                                                 / float(self.values.size))))
                col_val = float(wet_count)/colors
                dry_count = colors - wet_count
                colors_undersea = plt.cm.bwr(np.linspace(1., 0., wet_count))
                colors_land = plt.cm.terrain(np.linspace(0.25, 1., dry_count))
                colors = np.vstack((colors_undersea, colors_land))
                cmap = LinearSegmentedColormap.from_list('cut_terrain', colors)
                wlevels = np.linspace(vmin, 0.0, wet_count, endpoint=False)
                dlevels = np.linspace(0.0, vmax, dry_count)
                levels = np.hstack((wlevels, dlevels))
        else:
            cmap = plt.cm.get_cmap(cmap)
            levels = np.linspace(vmin, vmax, colors)
            col_val = 0.
        if vmax > 0:
            if norm is None:
                norm = FixPointNormalize(sealevel=0.0, vmax=vmax, vmin=vmin,
                                         col_val=col_val)
        return cmap, norm, levels, col_val

    def _add_boundary(self, btype):
        if self.has_attribute(btype):
            self.remove_attribute(btype)

        if btype == 'ocean_boundaries':
            self.add_attribute(
                btype,
                **{'values': None,
                   'index': None})

        elif btype in ['land_boundaries', 'inner_boundaries',
                       'inflow_boundaries']:
            self.add_attribute(
                btype,
                **{'values': None,
                   'index': None,
                   'ibtype': None})

        elif btype == 'outflow_boundaries':
            self.add_attribute(
                btype,
                **{'values': None,
                   'ibtype': None,
                   'index': None,
                   'barrier_height': None,
                   'supercritical_flow_coefficient': None})

        elif btype == 'weir_boundaries':
            self.add_attribute(
                btype,
                **{'values': None,
                   'ibtype': None,
                   'front_face_index': None,
                   'back_face_index': None,
                   'barrier_height': None,
                   'subcritical_flow_coefficient': None,
                   'supercritical_flow_coefficient': None})

        elif btype == 'culvert_boundaries':
            self.add_attribute(
                btype,
                **{'values': None,
                   'ibtype': None,
                   'front_face_index': None,
                   'back_face_index': None,
                   'barrier_height': None,
                   'subcritical_flow_coefficient': None,
                   'supercritical_flow_coefficient': None,
                   'cross_barrier_pipe_height': None,
                   'friction_factor': None,
                   'pipe_diameter': None})

    def _add_boundary_data(self, btype, bdata):

        mask = np.full((self.values.shape[0],), -99999, dtype=int)
        for i, boundary in enumerate(bdata):
            if btype in ['weir_boundaries', 'culvert_boundaries']:
                mask[boundary['front_face_indexes']] = i
                mask[boundary['back_face_indexes']] = -i
            else:
                mask[boundary['indexes']] = i

        boundary_collection = dict()
        if btype in 'ocean_boundaries':
            boundary_collection['index'] = list()
            for boundary in bdata:
                boundary_collection['index'].append(boundary['indexes'])

        if btype in ['land_boundaries', 'inner_boundaries',
                     'inflow_boundaries']:
            boundary_collection['index'] = list()
            boundary_collection['ibtype'] = list()
            for boundary in bdata:
                boundary_collection['index'].append(boundary['indexes'])
                boundary_collection['ibtype'].append(boundary['ibtype'])

        if btype == 'outflow_boundaries':
            boundary_collection['index'] = list()
            boundary_collection['ibtype'] = list()
            boundary_collection['barrier_height'] = list()
            boundary_collection['supercritical_flow_coefficient'] = list()
            for boundary in bdata:
                boundary_collection['index'].append(boundary['indexes'])
                boundary_collection['ibtype'].append(boundary['ibtype'])
                boundary_collection['supercritical_flow_coefficient'].append(
                    boundary['supercritical_flow_coefficients'])
                boundary_collection['barrier_height'].append(
                    boundary['barrier_heights'])

        if btype == 'weir_boundaries':
            boundary_collection['front_face_index'] = list()
            boundary_collection['back_face_index'] = list()
            boundary_collection['ibtype'] = list()
            boundary_collection['barrier_height'] = list()
            boundary_collection['supercritical_flow_coefficient'] = list()
            boundary_collection['subcritical_flow_coefficient'] = list()
            for boundary in bdata:
                boundary_collection['front_face_index'].append(
                    boundary['front_face_indexes'])
                boundary_collection['back_face_index'].append(
                    boundary['back_face_indexes'])
                boundary_collection['ibtype'].append(boundary['ibtype'])
                boundary_collection['supercritical_flow_coefficient'].append(
                    boundary['supercritical_flow_coefficients'])
                boundary_collection['barrier_height'].append(
                    boundary['barrier_heights'])
                boundary_collection['subcritical_flow_coefficient'].append(
                    boundary['subcritical_flow_coefficients'])

        if btype == 'culvert_boundaries':
            boundary_collection['front_face_index'] = list()
            boundary_collection['back_face_index'] = list()
            boundary_collection['ibtype'] = list()
            boundary_collection['barrier_height'] = list()
            boundary_collection['supercritical_flow_coefficient'] = list()
            boundary_collection['subcritical_flow_coefficient'] = list()
            boundary_collection['cross_barrier_pipe_heights'] = list()
            boundary_collection['friction_factors'] = list()
            boundary_collection['pipe_diameters'] = list()
            for boundary in bdata:
                boundary_collection['front_face_index'].append(
                    boundary['front_face_indexes'])
                boundary_collection['back_face_index'].append(
                    boundary['back_face_indexes'])
                boundary_collection['ibtype'].append(boundary['ibtype'])
                boundary_collection['supercritical_flow_coefficient'].append(
                    boundary['supercritical_flow_coefficients'])
                boundary_collection['barrier_height'].append(
                    boundary['barrier_heights'])
                boundary_collection['subcritical_flow_coefficient'].append(
                    boundary['subcritical_flow_coefficients'])

        self.set_attribute(btype, mask, **boundary_collection)

    def _get_boundary_collection(self, btype):
        if btype in ['ocean_boundaries', 'land_boundaries', 'inner_boundaries',
                     'inflow_boundaries', 'outflow_boundaries']:
            return self.get_attribute(btype)['index']
        else:
            bdata = self.get_attribute(btype)
            index_collection = list()
            for i in range(len(bdata['front_face_index'])):
                index_collection.append(list(
                    zip(bdata['front_face_index'][i],
                        bdata['back_face_index'][i])))
            return index_collection

    @staticmethod
    def parse_fort14(path):
        grd = dict()
        grd['x'] = list()
        grd['y'] = list()
        grd['values'] = list()
        grd['node_id'] = list()
        grd['elements'] = list()
        grd['element_id'] = list()
        grd['ocean_boundaries'] = list()
        grd['land_boundaries'] = list()
        grd['inner_boundaries'] = list()
        grd['inflow_boundaries'] = list()
        grd['outflow_boundaries'] = list()
        grd['weir_boundaries'] = list()
        grd['culvert_boundaries'] = list()
        with open(pathlib.Path(path), 'r') as f:
            grd['description'] = "{}".format(f.readline())
            NE, NP = map(int, f.readline().split())
            for i in range(NP):
                node_id, x, y, z = f.readline().split()
                grd['node_id'].append(int(node_id)-1)
                grd['x'].append(float(x))
                grd['y'].append(float(y))
                grd['values'].append(-float(z))
            for i in range(NE):
                line = f.readline().split()
                grd['element_id'].append(int(line[0])-1)
                grd['elements'].append([int(x)-1 for x in line[2:]])
            msg = 'ERROR: only supports triangular meshes.'
            assert all(len(element) == 3 for element in grd['elements']), msg
            grd['vertices'] = np.vstack([grd.pop('x'), grd.pop('y')]).T
            # Assume EOF if NOPE is empty.
            try:
                NOPE = int(f.readline().split()[0])
            except IndexError:
                return grd
            # For now, let NOPE=-1 mean a self closing mesh
            # reassigning NOPE to 0 until further implementation is applied.
            if NOPE == -1:
                NOPE = 0
            _NOPE = len([])
            f.readline()  # Number of total open ocean nodes. Not used.
            while _NOPE < NOPE:
                grd['ocean_boundaries'].append({'indexes': list()})
                NETA = int(f.readline().split()[0])
                _NETA = len([])
                while _NETA < NETA:
                    grd['ocean_boundaries'][_NOPE]['indexes'].append(
                                                int(f.readline().split()[0])-1)
                    _NETA += 1
                _NOPE += 1
            NBOU = int(f.readline().split()[0])
            _NBOU = len([])
            f.readline()
            while _NBOU < NBOU:
                NVELL, IBTYPE = map(int, f.readline().split()[:2])
                _NVELL = 0
                if IBTYPE in [0, 10, 20, 52]:  # ADDED IBTYPE 52 HERE ALTOUGH THAT MIGHT NOT BE TECHNICALLY CORRECT
                    grd['land_boundaries'].append(
                        {'ibtype': IBTYPE,
                         'indexes': list()})
                elif IBTYPE in [1, 11, 21]:
                    grd['inner_boundaries'].append(
                        {'ibtype': IBTYPE,
                         'indexes': list()})
                elif IBTYPE in [2, 12, 22, 102, 122]:
                    grd['inflow_boundaries'].append(
                        {'ibtype': IBTYPE,
                         'indexes': list()})
                elif IBTYPE in [3, 13, 23]:
                    grd['outflow_boundaries'].append(
                        {'ibtype': IBTYPE,
                         'indexes': list(),
                         'barrier_heights': list(),
                         'supercritical_flow_coefficients': list()})

                elif IBTYPE in [4, 24]:
                    grd['weir_boundaries'].append(
                        {'ibtype': IBTYPE,
                         'front_face_indexes': list(),
                         'back_face_indexes': list(),
                         'barrier_heights': list(),
                         'subcritical_flow_coefficients': list(),
                         'supercritical_flow_coefficients': list()})
                elif IBTYPE in [5, 25]:
                    grd['culvert_boundaries'].append(
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
                    if IBTYPE in [0, 10, 20, 52]:  # ADDED IBTYPE 52 HERE ALTOUGH THAT MIGHT NOT BE TECHNICALLY CORRECT
                        grd['land_boundaries'][-1][
                            'indexes'].append(int(line[0])-1)
                    elif IBTYPE in [1, 11, 21]:
                        grd['inner_boundaries'][-1][
                            'indexes'].append(int(line[0])-1)
                    elif IBTYPE in [2, 12, 22, 102, 122]:
                        grd['inflow_boundaries'][-1][
                            'indexes'].append(int(line[0])-1)
                    elif IBTYPE in [3, 13, 23]:
                        grd['outflow_boundaries'][-1][
                            'indexes'].append(int(line[0])-1)
                        grd['outflow_boundaries'][-1][
                            'barrier_heights'].append(float(line[1]))
                        grd['outflow_boundaries'][-1][
                            'supercritical_flow_coefficients'].append(
                                float(line[2]))
                    elif IBTYPE in [4, 24]:
                        grd['weir_boundaries'][-1][
                            'front_face_indexes'].append(int(line[0])-1)
                        grd['weir_boundaries'][-1][
                            'back_face_indexes'].append(int(line[1])-1)
                        grd['weir_boundaries'][-1][
                            'barrier_heights'].append(float(line[2]))
                        grd['weir_boundaries'][-1][
                            'subcritical_flow_coefficients'].append(
                                float(line[3]))
                        grd['weir_boundaries'][-1][
                            'supercritical_flow_coefficients'].append(
                                float(line[4]))
                    elif IBTYPE in [5, 25]:
                        grd['culvert_boundaries'][-1][
                            'front_face_indexes'].append(int(line[0])-1)
                        grd['culvert_boundaries'][-1][
                            'back_face_indexes'].append(int(line[1])-1)
                        grd['culvert_boundaries'][-1][
                            'barrier_heights'].append(float(line[2]))
                        grd['culvert_boundaries'][-1][
                            'subcritical_flow_coefficients'].append(
                                float(line[3]))
                        grd['culvert_boundaries'][-1][
                            'supercritical_flow_coefficients'].append(
                                float(line[4]))
                        grd['culvert_boundaries'][-1][
                            'friction_factors'].append(float(line[5]))
                        grd['culvert_boundaries'][-1][
                            'pipe_diameters'].append(float(line[6]))
                    else:
                        Exception("Duck-typing error. "
                                  + "This exception should be unreachable.")
                    _NVELL += 1
                _NBOU += 1
        return grd

    @staticmethod
    def parse_fort13(path):
        fort13 = dict()
        with open(pathlib.Path(path), 'r') as f:
            fort13['AGRID'] = f.readline()
            fort13['NumOfNodes'] = int(f.readline())
            num_of_attributes = int(f.readline().split()[0])
            i = 0
            while i < num_of_attributes:
                attribute_name = f.readline().strip()
                units = f.readline().strip()
                if units == '1':
                    units = 'unitless'
                f.readline()
                defaults = [float(x) for x in f.readline().split()]
                fort13[attribute_name] = dict()
                fort13[attribute_name]['units'] = units
                fort13[attribute_name]['default_values'] = defaults
                i += 1
            for i in range(num_of_attributes):
                attribute_name = f.readline().strip()
                num_of_nodes = int(f.readline())
                indexes = list()
                values = list()
                j = 0
                while j < num_of_nodes:
                    line = f.readline().split()
                    indexes.append(int(line.pop(0))-1)
                    line = [float(x) for x in line]
                    values.append(line)
                    j += 1
                fort13[attribute_name]['indexes'] = indexes
                values = np.asarray(values).flatten()
                # if values.shape[1] == 1:
                #     values = values.flatten()
                fort13[attribute_name]['values'] = values
        return fort13

    @staticmethod
    def _signed_polygon_area(vertices):
        # https://code.activestate.com/recipes/578047-area-of-polygon-using-shoelace-formula/
        n = len(vertices)  # of vertices
        area = 0.0
        for i in range(n):
            j = (i + 1) % n
            area += vertices[i][0] * vertices[j][1]
            area -= vertices[j][0] * vertices[i][1]
        return area / 2.0

    @property
    def description(self):
        return self._description

    @property
    def values(self):
        return self._values

    @property
    def node_id(self):
        return self._node_id

    @property
    def element_id(self):
        return self._element_id

    @property
    def ocean_boundaries(self):
        return self._ocean_boundaries

    @property
    def land_boundaries(self):
        return self._land_boundaries

    @property
    def inner_boundaries(self):
        return self._inner_boundaries

    @property
    def inflow_boundaries(self):
        return self._inflow_boundaries

    @property
    def outflow_boundaries(self):
        return self._outflow_boundaries

    @property
    def weir_boundaries(self):
        return self._weir_boundaries

    @property
    def culvert_boundaries(self):
        return self._culvert_boundaries

    @property
    def nodal_attribute_collection(self):
        return self._nodal_attribute_collection.copy()

    @property
    def primitive_weighting_in_continuity_equation(self):
        return self.get_attribute(
            "primitive_weighting_in_continuity_equation")

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

    @property
    def chezy_friction_coefficient_at_sea_floor(self):
        return self.get_attribute(
            "chezy_friction_coefficient_at_sea_floor")

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
    def fort14(self):
        f = f"{self.description}\n"
        f += f"{self.NE} {self.NP}\n"
        for i in range(self.NP):
            f += f"{i+1:d} "
            f += f"{self.x[i]:<.16E} "
            f += f"{self.y[i]:<.16E} "
            f += f"{-self.values[i]:<.16E}\n"
        for i in range(self.NE):
            f += f"{i+1:d} "
            f += f"{len(self.elements[i])} "
            f += f"{self.elements[i, 0]+1:d} "
            f += f"{self.elements[i, 1]+1:d} "
            f += f"{self.elements[i, 2]+1:d}\n"
        # ocean boundaries
        f += f"{len(self.ocean_boundaries):d} "
        f += "! total number of ocean boundaries\n"
        # count total number of ocean boundaries
        _sum = np.sum([len(boundary) for boundary in self.ocean_boundaries])
        f += f"{int(_sum):d} ! total number of ocean boundary nodes\n"
        # write ocean boundary indexes
        for i, boundary in enumerate(self.ocean_boundaries):
            f += f"{len(boundary):d}"
            f += f" ! number of nodes for ocean_boundary_{i}\n"
            for idx in boundary:
                f += f"{idx+1:d}\n"
        # count remaining boundaries
        f += "{:d}".format(
            len(self.land_boundaries) +
            len(self.inner_boundaries) +
            len(self.inflow_boundaries) +
            len(self.outflow_boundaries) +
            len(self.weir_boundaries) +
            len(self.culvert_boundaries))
        # count total remaining boundary simplices
        f += " ! total number of non-ocean boundaries\n"
        f += "{:d} ! Total number of non-ocean boundaries nodes\n".format(
            int(np.sum([len(x) for x in self.land_boundaries]) +
                np.sum([len(x) for x in self.inner_boundaries]) +
                np.sum([len(x) for x in self.inflow_boundaries]) +
                np.sum([len(x) for x in self.outflow_boundaries]) +
                np.sum([2*len(x) for x in self.weir_boundaries]) +
                np.sum([2*len(x) for x in self.culvert_boundaries])))
        # write land boundaries
        for i, boundary in enumerate(self.land_boundaries):
            ibtype = self.get_attribute('land_boundaries')['ibtype'][i]
            f += f"{len(boundary):d} "
            f += f"{ibtype:d} "
            f += "! number of nodes and ibtype for land_boundary_"
            f += f"{i}\n"
            for idx in boundary:
                f += f"{idx+1:d}\n"
        # write inner boundaries
        for i, boundary in enumerate(self.inner_boundaries):
            f += f"{len(boundary):d} "
            ibtype = self.get_attribute('inner_boundaries')['ibtype'][i]
            f += f"{ibtype:d} "
            f += "! number of nodes and ibtype for inner_boundary_"
            f += f"{i}\n"
            for idx in boundary:
                f += f"{idx+1:d}\n"
        # inflow boundaries
        for i, boundary in enumerate(self.inflow_boundaries):
            f += f"{len(boundary):d} "
            ibtype = self.get_attribute('inflow_boundaries')['ibtype'][i]
            f += f"{ibtype:d} "
            f += "! number of nodes and ibtype for inflow_boundary_"
            f += f"{i}\n"
            for idx in boundary:
                f += f"{idx+1:d}\n"
        # outflow boundaries
        for i, boundary in enumerate(self.outflow_boundaries):
            f += f"{len(boundary):d} "
            bdata = self.get_attribute('outflow_boundaries')
            ibtype = bdata['ibtype'][i]
            f += f"{ibtype:d} "
            f += "! number of nodes and ibtype for outflow_boundary_"
            f += f"{i}\n"
            for j, idx in enumerate(boundary):
                f += f"{idx+1:d} "
                f += f"{bdata['barrier_height'][i][j]:G} "
                f += f"{bdata['supercritical_flow_coefficient'][i][j]:G}\n"
        # weir boundaries
        for i, boundary in enumerate(self.weir_boundaries):
            f += f"{len(boundary):d} "
            bdata = self.get_attribute('weir_boundaries')
            ibtype = bdata['ibtype'][i]
            f += f"{ibtype:d} "
            f += "! number of nodes and ibtype for weir_boundary_"
            f += f"{i}\n"
            for j, idx in enumerate(boundary):
                f += f"{idx[0]+1:d} "
                f += f"{idx[1]+1:d} "
                f += f"{bdata['barrier_height'][i][j]:G} "
                f += f"{bdata['subcritical_flow_coefficient'][i][j]:G} "
                f += f"{bdata['supercritical_flow_coefficient'][i][j]:G}\n"
        # culvert boundaries
        for i, boundary in enumerate(self.culvert_boundaries):
            f += f"{len(boundary):d} "
            bdata = self.get_attribute('culvert_boundaries')
            ibtype = bdata['ibtype'][i]
            f += f"{ibtype:d} "
            f += "! number of nodes and ibtype for culvert_boundary_"
            f += f"{i}\n"
            for j, idx in enumerate(boundary):
                f += f"{idx[0]+1:d} "
                f += f"{idx[1]+1:d} "
                f += f"{bdata['barrier_height'][i][j]:G} "
                f += f"{bdata['subcritical_flow_coefficient'][i][j]:G} "
                f += f"{bdata['supercritical_flow_coefficient'][i][j]:G} "
                f += f"{bdata['cross_barrier_pipe_height'][i][j]:G} "
                f += f"{bdata['friction_factor'][i][j]:G} "
                f += f"{bdata['pipe_diameter'][i][j]:G}\n"
        # f += f"{}\n".format(self.SpatialReference.ExportToWkt())
        return f

    @property
    def fort13(self):
        f = "{}\n".format(self.description)
        f += "{}\n".format(len(self.node_id))
        f += "{}\n".format(len(self.get_nodal_attribute_names()))
        for name in self.get_nodal_attribute_names():
            attribute = self.get_nodal_attribute(name)
            f += f"{name}\n"
            f += "{}\n".format(attribute['units'])
            f += '{}\n'.format(len(attribute['defaults']))
            for n in attribute['defaults']:
                f += f'{n:<.16E} '
            f += '\n'
        for name in self.get_nodal_attribute_names():
            attribute = self.get_nodal_attribute(name)
            f += f"{name}\n"
            f += f"{len(attribute['non_default_indexes'])}\n"
            for idx in attribute['non_default_indexes']:
                f += f"{idx+1:d} "
                for val in attribute['values'][idx, :]:
                    f += f"{val:<.16E} "
                f += "\n"
        return f

    @property
    def NP(self):
        return self.vertices.shape[0]

    @property
    def NE(self):
        return self.elements.shape[0]

    @property
    def H0(self):
        try:
            return self.__H0
        except AttributeError:
            self.pad_lakes(H0=0.05)
            return self.__H0

    @mannings_n_at_sea_floor.setter
    def mannings_n_at_sea_floor(self, mannings_n_at_sea_floor):
        self.add_nodal_attribute('mannings_n_at_sea_floor', 'meters')
        self.set_nodal_attribute(
            'mannings_n_at_sea_floor', mannings_n_at_sea_floor)

    @property
    def _description(self):
        return self.__description

    @property
    def _nodal_attributes(self):
        return self.__nodal_attributes

    @property
    def _values(self):
        return self.get_attribute('elevation')['values']

    @property
    def _node_id(self):
        return self.__node_id

    @property
    def _element_id(self):
        return self.__element_id

    @property
    def _ocean_boundaries(self):
        return self._get_boundary_collection('ocean_boundaries')

    @property
    def _land_boundaries(self):
        return self._get_boundary_collection('land_boundaries')

    @property
    def _inner_boundaries(self):
        return self._get_boundary_collection('inner_boundaries')

    @property
    def _inflow_boundaries(self):
        return self._get_boundary_collection('inflow_boundaries')

    @property
    def _outflow_boundaries(self):
        return self._get_boundary_collection('outflow_boundaries')

    @property
    def _weir_boundaries(self):
        return self._get_boundary_collection('weir_boundaries')

    @property
    def _culvert_boundaries(self):
        return self._get_boundary_collection('culvert_boundaries')

    @property
    def _nodal_attribute_collection(self):
        try:
            return self.__nodal_attribute_collection
        except AttributeError:
            self.__nodal_attribute_collection = dict()
            return self.__nodal_attribute_collection

    @property
    def _fort13(self):
        return self.__fort13

    @_description.setter
    def _description(self, description):
        if description is None:
            description = uuid.uuid4().hex[:8]
        assert isinstance(description, str)
        msg = 'description field is limited to 32 characters.'
        description = description.strip('\n').strip()
        assert len(description) <= 32, msg
        self.__description = description

    @_node_id.setter
    def _node_id(self, node_id):
        if node_id is None:
            node_id = np.arange(self.vertices.shape[0])
        else:
            node_id = np.asarray(node_id).flatten()
            msg = "ERROR: Subgrids are yet not implemented."
            assert np.array_equal(
                node_id, np.arange(self.vertices.shape[0])), msg
        self.__node_id = node_id

    @_element_id.setter
    def _element_id(self, element_id):
        if element_id is None:
            element_id = np.arange(self.elements.shape[0])
        else:
            element_id = np.asarray(element_id).flatten()
            msg = "ERROR: Subgrids are yet not implemented."
            assert np.array_equal(
                element_id, np.arange(self.elements.shape[0])), msg
        self.__element_id = element_id

    @_values.setter
    def _values(self, values):
        values = np.array(values).flatten()
        assert values.size == self.vertices.shape[0]
        if self.has_attribute('elevation'):
            self.remove_attribute('elevation')
        self.add_attribute('elevation')
        self.set_attribute('elevation', values)

    @_ocean_boundaries.setter
    def _ocean_boundaries(self, ocean_boundaries):
        self._add_boundary('ocean_boundaries')
        ocean_boundaries = [] if ocean_boundaries is None else ocean_boundaries
        for i, ocean_boundary in enumerate(ocean_boundaries):
            assert 'indexes' in ocean_boundary.keys()
        self._add_boundary_data('ocean_boundaries', ocean_boundaries)

    @_land_boundaries.setter
    def _land_boundaries(self, land_boundaries):
        self._add_boundary('land_boundaries')
        land_boundaries = [] if land_boundaries is None else land_boundaries
        for i, land_boundary in enumerate(land_boundaries):
            assert 'indexes' in land_boundary.keys()
            assert 'ibtype' in land_boundary.keys()
            assert land_boundary['ibtype'] in [0, 10, 20, 52]  # ADDED IBTYPE 52 HERE ALTOUGH THAT MIGHT NOT BE TECHNICALLY CORRECT
        self._add_boundary_data('land_boundaries', land_boundaries)

    @_inner_boundaries.setter
    def _inner_boundaries(self, inner_boundaries):
        self._add_boundary('inner_boundaries')
        inner_boundaries = [] if inner_boundaries is None else inner_boundaries
        for i, inner_boundary in enumerate(inner_boundaries):
            assert 'indexes' in inner_boundary.keys()
            assert 'ibtype' in inner_boundary.keys()
            assert inner_boundary['ibtype'] in [1, 11, 21]
        self._add_boundary_data('inner_boundaries', inner_boundaries)

    @_inflow_boundaries.setter
    def _inflow_boundaries(self, inflow_boundaries):
        self._add_boundary('inflow_boundaries')
        inflow_boundaries = [] if inflow_boundaries is None \
            else inflow_boundaries
        for i, inflow_boundary in enumerate(inflow_boundaries):
            assert 'indexes' in inflow_boundary.keys()
            assert 'ibtype' in inflow_boundary.keys()
            assert inflow_boundary['ibtype'] in [2, 12, 22, 102, 122]
        self._add_boundary_data('inflow_boundaries', inflow_boundaries)

    @_outflow_boundaries.setter
    def _outflow_boundaries(self, outflow_boundaries):
        self._add_boundary('outflow_boundaries')
        outflow_boundaries = [] if outflow_boundaries is None \
            else outflow_boundaries
        for i, outflow_boundary in enumerate(outflow_boundaries):
            assert 'indexes' in outflow_boundary.keys()
            assert 'ibtype' in outflow_boundary.keys()
            assert 'barrier_heights' in outflow_boundary.keys()
            assert 'supercritical_flow_coefficients' in outflow_boundary.keys()
            assert outflow_boundary['ibtype'] in [3, 13, 23]
        self._add_boundary_data('outflow_boundaries', outflow_boundaries)

    @_weir_boundaries.setter
    def _weir_boundaries(self, weir_boundaries):
        self._add_boundary('weir_boundaries')
        weir_boundaries = [] if weir_boundaries is None \
            else weir_boundaries
        for i, weir_boundary in enumerate(weir_boundaries):
            assert 'front_face_indexes' in weir_boundary.keys()
            assert 'back_face_indexes' in weir_boundary.keys()
            assert 'ibtype' in weir_boundary.keys()
            assert 'barrier_heights' in weir_boundary.keys()
            assert 'supercritical_flow_coefficients' in weir_boundary.keys()
            assert 'subcritical_flow_coefficients' in weir_boundary.keys()
            assert weir_boundary['ibtype'] in [4, 24]
        self._add_boundary_data('weir_boundaries', weir_boundaries)

    @_culvert_boundaries.setter
    def _culvert_boundaries(self, culvert_boundaries):
        self._add_boundary('culvert_boundaries')
        culvert_boundaries = [] if culvert_boundaries is None \
            else culvert_boundaries
        for i, culvert_boundary in enumerate(culvert_boundaries):
            assert 'front_face_indexes' in culvert_boundary.keys()
            assert 'back_face_indexes' in culvert_boundary.keys()
            assert 'ibtype' in culvert_boundary.keys()
            assert 'barrier_heights' in culvert_boundary.keys()
            assert 'supercritical_flow_coefficients' in culvert_boundary.keys()
            assert 'subcritical_flow_coefficients' in culvert_boundary.keys()
            assert 'cross_barrier_pipe_heights' in culvert_boundary.keys()
            assert 'friction_factors' in culvert_boundary.keys()
            assert 'pipe_diameters' in culvert_boundary.keys()
            assert culvert_boundary['ibtype'] in [5, 25]
        self._add_boundary_data('culvert_boundaries', culvert_boundaries)

    @_fort13.setter
    def _fort13(self, fort13):
        if fort13 is not None:
            self.import_nodal_attributes(fort13)
        self.__fort13 = fort13


#     def add_ocean_boundary(self, key, indexes):
#         if key in self.__ocean_boundaries.keys():
#             raise AttributeError(
#                 'Cannot add ocean boundary with key {}:'.format(key)
#                 + ' key already exists.')
#         else:
#             self.__ocean_boundaries[key] = dict()
#         ocean_boundary = np.full(self.values.shape, False)
#         ocean_boundary[np.where(np.in1d(self.node_id, indexes))] = True
#         self.__ocean_boundaries[key]['__bool'] = ocean_boundary

#     def add_land_boundary(self, key, indexes, ibtype=20):
#         assert ibtype in [0, 10, 20]
#         if key in self.__land_boundaries.keys():
#             raise AttributeError(
#                 'Cannot add land boundary with key {}:'.format(key)
#                 + ' key already exists.')
#         else:
#             self.__land_boundaries[key] = dict()
#         land_boundary = np.full(self.values.shape, False)
#         land_boundary[np.where(np.in1d(self.node_id, indexes))] = True
#         self.__land_boundaries[key]['__bool'] = land_boundary
#         self.__land_boundaries[key]['ibtype'] = ibtype

#     def add_inner_boundary(self, key, indexes, ibtype=21):
#         assert ibtype in [1, 11, 21]
#         if key in self.__inner_boundaries.keys():
#             raise AttributeError(
#                 'Cannot add inner boundary with key {}:'.format(key)
#                 + ' key already exists.')
#         else:
#             self.__inner_boundaries[key] = dict()
#         inner_boundary = np.full(self.values.shape, False)
#         inner_boundary[np.where(np.in1d(self.node_id, indexes))] = True
#         self.__inner_boundaries[key]['__bool'] = inner_boundary
#         self.__inner_boundaries[key]['ibtype'] = ibtype

#     def add_inflow_boundary(self, key, indexes, ibtype=22):
#         assert ibtype in [2, 12, 22, 102, 122]
#         if key in self.__inflow_boundaries.keys():
#             raise AttributeError(
#                 'Cannot add inflow boundary with key {}:'.format(key)
#                 + ' key already exists.')
#         else:
#             self.__inflow_boundaries[key] = dict()
#         inflow_boundary = np.full(self.values.shape, False)
#         inflow_boundary[np.where(np.in1d(self.node_id, indexes))] = True
#         self.__inflow_boundaries[key]['__bool'] = inflow_boundary
#         self.__inflow_boundaries[key]['ibtype'] = ibtype

#     def add_outflow_boundary(
#         self,
#         key,
#         indexes,
#         barrier_heights,
#         supercritical_flow_coefficients,
#         ibtype=23
#     ):
#         assert ibtype in [3, 13, 23]
#         indexes = np.asarray(indexes)
#         if key in self.__outflow_boundaries.keys():
#             raise AttributeError(
#                 'Cannot add outflow boundary with key {}:'.format(key)
#                 + ' key already exists.')
#         else:
#             self.__outflow_boundaries[key] = dict()
#         outflow_boundary = np.full(self.values.shape, False)
#         outflow_boundary[np.where(np.in1d(self.node_id, indexes))] = True
#         barrier_heights = np.asarray(barrier_heights)
#         assert barrier_heights.size == indexes.size
#         supercritical_flow_coefficient \
#             = np.asarray(supercritical_flow_coefficients)
#         assert supercritical_flow_coefficient.size == indexes.size
#         self.__outflow_boundaries[key]['__bool'] = outflow_boundary
#         self.__outflow_boundaries[key]['ibtype'] = ibtype
#         self.__outflow_boundaries[key]['barrier_heights'] = barrier_heights
#         self.__outflow_boundaries[key]['supercritical_flow_coefficients'] \
#             = supercritical_flow_coefficient

#     def add_weir_boundary(
#         self,
#         key,
#         front_face_indexes,
#         back_face_indexes,
#         barrier_heights,
#         subcritical_flow_coefficients,
#         supercritical_flow_coefficients,
#         ibtype=24
#     ):
#         assert ibtype in [4, 24]
#         front_face_indexes = np.asarray(front_face_indexes).flatten()
#         back_face_indexes = np.asarray(back_face_indexes).flatten()
#         assert front_face_indexes.size == back_face_indexes.size
#         if key in self.__weir_boundaries.keys():
#             raise AttributeError(
#                 'Cannot add weir boundary with key {}:'.format(key)
#                 + ' key already exists.')
#         else:
#             self.__weir_boundaries[key] = dict()
#         weir_boundary_front_face = np.full(self.values.shape, False)
#         weir_boundary_front_face[np.where(np.in1d(self.node_id,
#                                  front_face_indexes))] = True
#         weir_boundary_back_face = np.full(self.values.shape, False)
#         weir_boundary_back_face[np.where(np.in1d(self.node_id,
#                                 back_face_indexes))] = True
#         barrier_heights = np.asarray(barrier_heights)
#         assert barrier_heights.size == front_face_indexes.size
#         subcritical_flow_coefficient \
#             = np.asarray(subcritical_flow_coefficients)
#         assert subcritical_flow_coefficient.size == front_face_indexes.size
#         supercritical_flow_coefficient \
#             = np.asarray(supercritical_flow_coefficients)
#         assert supercritical_flow_coefficient.size == front_face_indexes.size
#         self.__weir_boundaries[key]['__bool_front_face'] = \
#             weir_boundary_front_face
#         self.__weir_boundaries[key]['__bool_back_face'] = \
#             weir_boundary_back_face
#         self.__weir_boundaries[key]['ibtype'] = ibtype
#         self.__weir_boundaries[key]['barrier_heights'] = barrier_heights
#         self.__weir_boundaries[key]['subcritical_flow_coefficients'] \
#             = subcritical_flow_coefficient
#         self.__weir_boundaries[key]['supercritical_flow_coefficients'] \
#             = supercritical_flow_coefficient

#     def add_culvert_boundary(
#         self,
#         key,
#         front_face_indexes,
#         back_face_indexes,
#         barrier_heights, subcritical_flow_coefficients,
#         supercritical_flow_coefficients,
#         cross_barrier_pipe_heights, friction_factors,
#         pipe_diameters,
#         ibtype=25
#     ):
#         assert ibtype in [5, 25]
#         front_face_indexes = np.asarray(front_face_indexes).flatten()
#         back_face_indexes = np.asarray(back_face_indexes).flatten()
#         assert front_face_indexes.size == back_face_indexes.size
#         if key in self.__culvert_boundaries.keys():
#             raise AttributeError(
#                 'Cannot add culvert boundary with key {}:'.format(key)
#                 + ' key already exists.')
#         else:
#             self.__culvert_boundaries[key] = dict()
#         culvert_boundary_front_face = np.full(self.values.shape, False)
#         culvert_boundary_front_face[np.where(np.in1d(self.node_id,
#                                     front_face_indexes))] = True
#         culvert_boundary_back_face = np.full(self.values.shape, False)
#         culvert_boundary_back_face[np.where(np.in1d(self.node_id,
#                                    back_face_indexes))] = True
#         barrier_heights = np.asarray(barrier_heights)
#         assert barrier_heights.size == front_face_indexes.size
#         subcritical_flow_coefficient \
#             = np.asarray(subcritical_flow_coefficients)
#         assert subcritical_flow_coefficient.size == front_face_indexes.size
#         supercritical_flow_coefficient \
#             = np.asarray(supercritical_flow_coefficients)
#         assert supercritical_flow_coefficient.size == front_face_indexes.size
#         cross_barrier_pipe_heights = np.asarray(cross_barrier_pipe_heights)
#         assert cross_barrier_pipe_heights.size == front_face_indexes.size
#         friction_factors = np.asarray(friction_factors)
#         assert friction_factors.size == front_face_indexes.size
#         pipe_diameters = np.asarray(pipe_diameters)
#         assert pipe_diameters.size == front_face_indexes.size
#         self.__culvert_boundaries[key]['__bool_front_face'] = \
#             culvert_boundary_front_face
#         self.__culvert_boundaries[key]['__bool_back_face'] = \
#             culvert_boundary_back_face
#         self.__culvert_boundaries[key]['ibtype'] = ibtype
#         self.__culvert_boundaries[key]['barrier_heights'] = barrier_heights
#         self.__culvert_boundaries[key]['subcritical_flow_coefficients'] \
#             = subcritical_flow_coefficient
#         self.__culvert_boundaries[key]['supercritical_flow_coefficients'] \
#             = supercritical_flow_coefficient
#         self.__culvert_boundaries[key]['cross_barrier_pipe_heights'] \
#             = cross_barrier_pipe_heights
#         self.__culvert_boundaries[key]['friction_factors'] = friction_factors
#         self.__culvert_boundaries[key]['pipe_diameters'] = pipe_diameters

#     def remove_ocean_boundary(self, key):
#         self.__remove_boundary(key, self.__ocean_boundaries, "ocean")

#     def remove_land_boundary(self, key):
#         self.__remove_boundary(key, self.__land_boundaries, "land")

#     def remove_inner_boundary(self, key):
#         self.__remove_boundary(key, self.__inner_boundaries, "inner")

#     def remove_inflow_boundary(self, key):
#         self.__remove_boundary(key, self.__inflow_boundaries, "inflow")

#     def remove_outflow_boundary(self, key):
#         self.__remove_boundary(key, self.__outflow_boundaries, "outflow")

#     def remove_weir_boundary(self, key):
#         self.__remove_boundary(key, self.__weir_boundaries, "weir")

#     def remove_culvert_boundary(self, key):
#         self.__remove_boundary(key, self.__culvert_boundaries, "culvert")


    # def _get_boundary_collection(self, boundary_type):
    #     data = self.boundary_collection[boundary_type]
    #     if len(data.keys()) == 0:
    #         msg = 'No boundary information present for boundary type '
    #         msg += f'{boundary_type}.'
    #         raise AttributeError(msg)
    #     return data


        ##
        # values = self.get_attribute(boundary_type)['values']
        # boundary = np.ma.masked_equal(values, -99999)
        # boundary_collection = list()
        # for boundary_id in np.unique(boundary):
        #     if not isinstance(boundary_id, np.ma.core.MaskedConstant):
        #         # handles most boundary types
        #         if values.ndim == 1:
        #             idxs = np.where(boundary == boundary_id)[0]
        #             boundary_collection.append(
        #                 self._get_ordered_indexes(idxs))
        #         # handles boundary duals (weir, culvert)
        #         elif values.ndim == 2:
        #             _b = list()
        #             for i in range(values.ndim):
        #                 idxs = np.where(boundary[:, i] == boundary_id)[0]
        #                 _b.append(self._get_ordered_indexes(idxs))
        #             boundary_collection.append(_b)
        #         else:
        #             raise Exception('Boundary must be a simplex or dual.')
        # return boundary_collection

    # def _get_ordered_indexes(self, indexes):
    #     triangulation = self.triangulation
    #     boundary_edges = list()
    #     idxs = np.vstack(list(np.where(triangulation.neighbors == -1))).T
    #     for i, j in idxs:
    #         boundary_edges.append((triangulation.triangles[i, j],
    #                                triangulation.triangles[i, (j+1) % 3]))
    #     boundary_edges = np.asarray(boundary_edges)
    #     idxs = np.where(np.in1d(boundary_edges[:, 0], indexes))[0]
    #     boundary_edges = boundary_edges[idxs, :]
    #     _, unique = np.unique(boundary_edges.flatten(), return_counts=True)
    #     end_points = np.where(unique == 1)[0]
    #     if len(end_points) == 0:
    #         boundary_edges = boundary_edges.tolist()
    #         ordered_edges = [boundary_edges.pop(-1)]
    #     else:
    #         end_points = end_points % boundary_edges.shape[0]
    #         boundary_edges = boundary_edges.tolist()
    #         ordered_edges = [boundary_edges.pop(np.max(end_points))]
    #     while len(boundary_edges) > 0:
    #         current_length = len(boundary_edges)  # for avoiding infinite loop
    #         try:
    #             idx = np.where(
    #                 np.asarray(
    #                     boundary_edges)[:, 0] == ordered_edges[-1][1])[0][0]
    #             ordered_edges.append(boundary_edges.pop(idx))
    #         except IndexError:
    #             idx = np.where(
    #                 np.asarray(
    #                     boundary_edges)[:, 1] == ordered_edges[0][0])[0][0]
    #             ordered_edges.insert(0, boundary_edges.pop(idx))
    #         # raise on infinite loop
    #         msg = 'ERROR: Went into infinite loop when reading boundary data.'
    #         assert len(boundary_edges) == current_length - 1, msg
    #     return np.asarray(ordered_edges)[:, 0]