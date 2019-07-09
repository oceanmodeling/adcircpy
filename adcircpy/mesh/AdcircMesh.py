from copy import deepcopy
from collections import OrderedDict
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from osgeo import osr
from datetime import datetime
from collections import defaultdict
from itertools import permutations
from matplotlib.tri import Triangulation
from adcircpy.lib._FixPointNormalize import _FixPointNormalize
from adcircpy.mesh import UnstructuredMesh


class AdcircMesh(UnstructuredMesh):

    def __init__(
        self,
        vertices,
        elements,
        values,
        SpatialReference=None,
        node_id=None,
        element_id=None
    ):
        super(AdcircMesh, self).__init__(
            vertices, elements, values, SpatialReference)
        self.node_id = node_id
        self.element_id = element_id
        self.__ocean_boundaries = OrderedDict()
        self.__land_boundaries = OrderedDict()
        self.__inner_boundaries = OrderedDict()
        self.__inflow_boundaries = OrderedDict()
        self.__outflow_boundaries = OrderedDict()
        self.__weir_boundaries = OrderedDict()
        self.__culvert_boundaries = OrderedDict()
        self.__nodal_attributes = dict()

    def add_nodal_attribute(
        self,
        attribute,
        values,
        units=None,
        coldstart=False,
        hotstart=False
    ):
        self.add_attribute(attribute)
        self.set_attribute(attribute, values)
        if attribute in self.__nodal_attributes.keys():
            raise AttributeError(
                'Cannot add nodal attribute with name {}:'.format(attribute)
                + ' attribute already exists.')
        else:
            self.__nodal_attributes[attribute] = dict()
        assert isinstance(coldstart, bool)
        assert isinstance(hotstart, bool)
        self.__nodal_attributes[attribute] = {'units': units,
                                              'coldstart': coldstart,
                                              'hotstart': hotstart}

    def add_ocean_boundary(self, key, indexes):
        if key in self.__ocean_boundaries.keys():
            raise AttributeError(
                'Cannot add ocean boundary with key {}:'.format(key)
                + ' key already exists.')
        else:
            self.__ocean_boundaries[key] = dict()
        ocean_boundary = np.full(self.values.shape, False)
        ocean_boundary[np.where(np.in1d(self.node_id, indexes))] = True
        self.__ocean_boundaries[key]['__bool'] = ocean_boundary

    def add_land_boundary(self, key, indexes, ibtype=20):
        assert ibtype in [0, 10, 20]
        if key in self.__land_boundaries.keys():
            raise AttributeError(
                'Cannot add land boundary with key {}:'.format(key)
                + ' key already exists.')
        else:
            self.__land_boundaries[key] = dict()
        land_boundary = np.full(self.values.shape, False)
        land_boundary[np.where(np.in1d(self.node_id, indexes))] = True
        self.__land_boundaries[key]['__bool'] = land_boundary
        self.__land_boundaries[key]['ibtype'] = ibtype

    def add_inner_boundary(self, key, indexes, ibtype=21):
        assert ibtype in [1, 11, 21]
        if key in self.__inner_boundaries.keys():
            raise AttributeError(
                'Cannot add inner boundary with key {}:'.format(key)
                + ' key already exists.')
        else:
            self.__inner_boundaries[key] = dict()
        inner_boundary = np.full(self.values.shape, False)
        inner_boundary[np.where(np.in1d(self.node_id, indexes))] = True
        self.__inner_boundaries[key]['__bool'] = inner_boundary
        self.__inner_boundaries[key]['ibtype'] = ibtype

    def add_inflow_boundary(self, key, indexes, ibtype=22):
        assert ibtype in [2, 12, 22, 102, 122]
        if key in self.__inflow_boundaries.keys():
            raise AttributeError(
                'Cannot add inflow boundary with key {}:'.format(key)
                + ' key already exists.')
        else:
            self.__inflow_boundaries[key] = dict()
        inflow_boundary = np.full(self.values.shape, False)
        inflow_boundary[np.where(np.in1d(self.node_id, indexes))] = True
        self.__inflow_boundaries[key]['__bool'] = inflow_boundary
        self.__inflow_boundaries[key]['ibtype'] = ibtype

    def add_outflow_boundary(
        self,
        key,
        indexes,
        barrier_heights,
        supercritical_flow_coefficients,
        ibtype=23
    ):
        assert ibtype in [3, 13, 23]
        indexes = np.asarray(indexes)
        if key in self.__outflow_boundaries.keys():
            raise AttributeError(
                'Cannot add outflow boundary with key {}:'.format(key)
                + ' key already exists.')
        else:
            self.__outflow_boundaries[key] = dict()
        outflow_boundary = np.full(self.values.shape, False)
        outflow_boundary[np.where(np.in1d(self.node_id, indexes))] = True
        barrier_heights = np.asarray(barrier_heights)
        assert barrier_heights.size == indexes.size
        supercritical_flow_coefficient \
            = np.asarray(supercritical_flow_coefficients)
        assert supercritical_flow_coefficient.size == indexes.size
        self.__outflow_boundaries[key]['__bool'] = outflow_boundary
        self.__outflow_boundaries[key]['ibtype'] = ibtype
        self.__outflow_boundaries[key]['barrier_heights'] = barrier_heights
        self.__outflow_boundaries[key]['supercritical_flow_coefficients'] \
            = supercritical_flow_coefficient

    def add_weir_boundary(
        self,
        key,
        front_face_indexes,
        back_face_indexes,
        barrier_heights,
        subcritical_flow_coefficients,
        supercritical_flow_coefficients,
        ibtype=24
    ):
        assert ibtype in [4, 24]
        front_face_indexes = np.asarray(front_face_indexes).flatten()
        back_face_indexes = np.asarray(back_face_indexes).flatten()
        assert front_face_indexes.size == back_face_indexes.size
        if key in self.__weir_boundaries.keys():
            raise AttributeError(
                'Cannot add weir boundary with key {}:'.format(key)
                + ' key already exists.')
        else:
            self.__weir_boundaries[key] = dict()
        weir_boundary_front_face = np.full(self.values.shape, False)
        weir_boundary_front_face[np.where(np.in1d(self.node_id,
                                 front_face_indexes))] = True
        weir_boundary_back_face = np.full(self.values.shape, False)
        weir_boundary_back_face[np.where(np.in1d(self.node_id,
                                back_face_indexes))] = True
        barrier_heights = np.asarray(barrier_heights)
        assert barrier_heights.size == front_face_indexes.size
        subcritical_flow_coefficient \
            = np.asarray(subcritical_flow_coefficients)
        assert subcritical_flow_coefficient.size == front_face_indexes.size
        supercritical_flow_coefficient \
            = np.asarray(supercritical_flow_coefficients)
        assert supercritical_flow_coefficient.size == front_face_indexes.size
        self.__weir_boundaries[key]['__bool_front_face'] = \
            weir_boundary_front_face
        self.__weir_boundaries[key]['__bool_back_face'] = \
            weir_boundary_back_face
        self.__weir_boundaries[key]['ibtype'] = ibtype
        self.__weir_boundaries[key]['barrier_heights'] = barrier_heights
        self.__weir_boundaries[key]['subcritical_flow_coefficients'] \
            = subcritical_flow_coefficient
        self.__weir_boundaries[key]['supercritical_flow_coefficients'] \
            = supercritical_flow_coefficient

    def add_culvert_boundary(
        self,
        key,
        front_face_indexes,
        back_face_indexes,
        barrier_heights, subcritical_flow_coefficients,
        supercritical_flow_coefficients,
        cross_barrier_pipe_heights, friction_factors,
        pipe_diameters,
        ibtype=25
    ):
        assert ibtype in [5, 25]
        front_face_indexes = np.asarray(front_face_indexes).flatten()
        back_face_indexes = np.asarray(back_face_indexes).flatten()
        assert front_face_indexes.size == back_face_indexes.size
        if key in self.__culvert_boundaries.keys():
            raise AttributeError(
                'Cannot add culvert boundary with key {}:'.format(key)
                + ' key already exists.')
        else:
            self.__culvert_boundaries[key] = dict()
        culvert_boundary_front_face = np.full(self.values.shape, False)
        culvert_boundary_front_face[np.where(np.in1d(self.node_id,
                                    front_face_indexes))] = True
        culvert_boundary_back_face = np.full(self.values.shape, False)
        culvert_boundary_back_face[np.where(np.in1d(self.node_id,
                                   back_face_indexes))] = True
        barrier_heights = np.asarray(barrier_heights)
        assert barrier_heights.size == front_face_indexes.size
        subcritical_flow_coefficient \
            = np.asarray(subcritical_flow_coefficients)
        assert subcritical_flow_coefficient.size == front_face_indexes.size
        supercritical_flow_coefficient \
            = np.asarray(supercritical_flow_coefficients)
        assert supercritical_flow_coefficient.size == front_face_indexes.size
        cross_barrier_pipe_heights = np.asarray(cross_barrier_pipe_heights)
        assert cross_barrier_pipe_heights.size == front_face_indexes.size
        friction_factors = np.asarray(friction_factors)
        assert friction_factors.size == front_face_indexes.size
        pipe_diameters = np.asarray(pipe_diameters)
        assert pipe_diameters.size == front_face_indexes.size
        self.__culvert_boundaries[key]['__bool_front_face'] = \
            culvert_boundary_front_face
        self.__culvert_boundaries[key]['__bool_back_face'] = \
            culvert_boundary_back_face
        self.__culvert_boundaries[key]['ibtype'] = ibtype
        self.__culvert_boundaries[key]['barrier_heights'] = barrier_heights
        self.__culvert_boundaries[key]['subcritical_flow_coefficients'] \
            = subcritical_flow_coefficient
        self.__culvert_boundaries[key]['supercritical_flow_coefficients'] \
            = supercritical_flow_coefficient
        self.__culvert_boundaries[key]['cross_barrier_pipe_heights'] \
            = cross_barrier_pipe_heights
        self.__culvert_boundaries[key]['friction_factors'] = friction_factors
        self.__culvert_boundaries[key]['pipe_diameters'] = pipe_diameters

    def import_nodal_attributes(self, fort13):
        fort13 = self.parse_fort13(fort13)
        for attribute, data in fort13.items():
            values = list()
            for i in range(len(data['default_values'])):
                _values = np.full((self.values.size,), np.nan)
                idxs = np.where(np.in1d(data['indexes'], self.node_id))
                for j, idx in enumerate(idxs):
                    _values[idx] = data['values'][i][j]
                _values[np.where(np.isnan(_values))] = \
                    data['default_values'][i]
                values.append(_values)
            values = np.asarray(values)
            if values.size//self.vertices.shape[0] == 1:
                values = values.flatten()
            else:
                values = values.reshape((self.vertices.shape[0],
                                         values.size//self.vertices.shape[0]))
            self.add_nodal_attribute(attribute, values, units=data['units'])

    def get_nodal_attribute_names(self):
        return list(self.__nodal_attributes.keys())

    def get_nodal_attribute(self, attribute):
        assert attribute in self.__nodal_attributes.keys(), \
            "Cannot get attribute '{}'".format(attribute) \
            + ": attribute does not exist."
        return {'values': self.get_attribute(attribute),
                **self.__nodal_attributes[attribute]}

    def get_coldstart_attributes(self):
        coldstart_attributes = dict()
        for attribute in self.get_nodal_attribute_names():
            attr = self.get_nodal_attribute(attribute)
            if attr['coldstart']:
                coldstart_attributes[attribute] = attr
        return coldstart_attributes

    def get_hotstart_attributes(self):
        hotstart_attributes = dict()
        for attribute in self.get_nodal_attribute_names():
            attr = self.get_nodal_attribute(attribute)
            if attr['coldstart']:
                hotstart_attributes[attribute] = attr
        return hotstart_attributes

    def set_nodal_attribute_coldstart_state(self, attribute, state):
        assert isinstance(state, bool)
        self.__nodal_attributes[attribute]['coldstart'] = state

    def set_nodal_attribute_hotstart_state(self, attribute, state):
        assert isinstance(state, bool)
        self.__nodal_attributes[attribute]['hotstart'] = state

    def set_nodal_attribute_state(
        self,
        attribute,
        coldstart_state,
        hotstart_state
    ):
        self.set_nodal_attribute_coldstart_state(attribute, coldstart_state)
        self.set_nodal_attribute_hotstart_state(attribute, hotstart_state)

    def tau0_gen(
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
        """
        msg = "Cannot compute TAU0 with nan depth values."
        assert not np.any(np.isnan(self.values)), msg
        msg = "Cannot compute TAU0 with no SpatialReference set."
        assert isinstance(self.SpatialReference, osr.SpatialReference), msg
        x = self.get_x(3395)
        y = self.get_y(3395)
        tri = Triangulation(x, y, self.elements)
        neighbors = defaultdict(set)
        for simplex in tri.triangles:
            for i, j in permutations(simplex, 2):
                neighbors[i].add(j)
        points = [tuple(p) for p in np.vstack([x, y]).T]
        values = np.full(self.values.shape, default_value)
        for k, v in neighbors.items():
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
        self.add_nodal_attribute(
            'primitive_weighting_in_continuity_equation', values, 'unitless',
            True, True)

    def has_primitive_weighting(self, runtype=None):
        if runtype is None:
            runtype = ['coldstart', 'hotstart']
        else:
            list(runtype)
        for _ in runtype:
            if runtype == 'coldstart':
                if 'primitive_weighting_in_continuity_equation' in \
                        self.get_coldstart_attributes():
                    return True
                else:
                    continue
            else:
                if 'primitive_weighting_in_continuity_equation' in \
                        self.get_hotstart_attributes():
                    return True
                else:
                    continue
        return False

    def remove_ocean_boundary(self, key):
        self.__remove_boundary(key, self.__ocean_boundaries, "ocean")

    def remove_land_boundary(self, key):
        self.__remove_boundary(key, self.__land_boundaries, "land")

    def remove_inner_boundary(self, key):
        self.__remove_boundary(key, self.__inner_boundaries, "inner")

    def remove_inflow_boundary(self, key):
        self.__remove_boundary(key, self.__inflow_boundaries, "inflow")

    def remove_outflow_boundary(self, key):
        self.__remove_boundary(key, self.__outflow_boundaries, "outflow")

    def remove_weir_boundary(self, key):
        self.__remove_boundary(key, self.__weir_boundaries, "weir")

    def remove_culvert_boundary(self, key):
        self.__remove_boundary(key, self.__culvert_boundaries, "culvert")

    def write_fort14(self, path=None):
        fort14 = "{}\n".format(self.description)
        fort14 += "{}  {}\n".format(self.num_elements, self.num_nodes)
        for i in range(self.num_nodes):
            fort14 += "{:d} ".format(self.node_id[i])
            fort14 += "{:<.16E} ".format(self.x[i])
            fort14 += " {:<.16E} ".format(self.y[i])
            fort14 += "{:<.16E}\n".format(-self.values[i])
        for i in range(self.num_elements):
            fort14 += "{:d} ".format(self.element_id[i])
            fort14 += "{:d} ".format(3)
            fort14 += "{:d} ".format(self.elements[i, 0]+1)
            fort14 += "{:d} ".format(self.elements[i, 1]+1)
            fort14 += "{:d}\n".format(self.elements[i, 2]+1)
        fort14 += "{:d} ! total number of ocean boundaries\n".format(
            len(self.ocean_boundaries.keys()))
        fort14 += "{:d} ! total number of ocean boundary nodes\n".format(
            len(self.ocean_boundary))
        for key, indexes in self.ocean_boundaries.items():
            fort14 += "{:d}".format(len(indexes))
            fort14 += " ! number of nodes for ocean_boundary_"
            fort14 += "{}\n".format(key)
            for idx in indexes:
                fort14 += "{:d}\n".format(idx+1)
        fort14 += "{:d}".format(
            len(self.land_boundaries.keys()) +
            len(self.inner_boundaries.keys()) +
            len(self.inflow_boundaries.keys()) +
            len(self.outflow_boundaries.keys()) +
            len(self.weir_boundaries.keys()) +
            len(self.culvert_boundaries.keys()))
        fort14 += " ! total number of non-ocean boundaries\n"
        fort14 += "{:d}".format(
            len(self.land_boundary) +
            len(self.inner_boundary) +
            len(self.inflow_boundary) +
            len(self.outflow_boundary) +
            len(self.weir_boundary) +
            len(self.culvert_boundary))
        fort14 += "! total number of non-ocean boundary nodes\n"
        for key, _ in self.land_boundaries.items():
            fort14 += "{:d} ".format(len(_['indexes']))
            fort14 += "{:d} ".format(_['ibtype'])
            fort14 += "! number of nodes and ibtype for land_boundary_"
            fort14 += "{}\n".format(key)
            for idx in _['indexes']:
                fort14 += "{:d}\n".format(idx+1)
        for key, _ in self.inner_boundaries.items():
            fort14 += "{:d} ".format(len(_['indexes']))
            fort14 += "{:d} ".format(_['ibtype'])
            fort14 += "! number of nodes and ibtype for inner_boundary_"
            fort14 += "{}\n".format(key)
            for idx in _['indexes']:
                fort14 += "{:d}\n".format(idx+1)
        for key, _ in self.inflow_boundaries.items():
            fort14 += "{:d} ".format(len(_['indexes']))
            fort14 += "{:d} ".format(_['ibtype'])
            fort14 += "! number of nodes and ibtype for inflow_boundary_"
            fort14 += "{}\n".format(key)
            for idx in _['indexes']:
                fort14 += "{:d}\n".format(idx+1)
        for key, _ in self.outflow_boundaries.items():
            fort14 += "{:d} ".format(len(_['indexes']))
            fort14 += "{:d} ".format(_['ibtype'])
            fort14 += "! number of nodes and ibtype for outflow_boundary_"
            fort14 += "{}\n".format(key)
            for i in range(len(_['indexes'])):
                fort14 += "{:d} ".format(_['indexes'][i]+1)
                fort14 += "{:<.16E} ".format(_["barrier_heights"][i])
                fort14 += "{:<.16E} ".format(
                        _["subcritical_flow_coefficients"][i])
                fort14 += "\n"
        for key, _ in self.weir_boundaries.items():
            fort14 += "{:d} ".format(len(_['front_face_indexes']))
            fort14 += "{:d} ".format(_['ibtype'])
            fort14 += "! number of nodes and ibtype for weir_boundary_"
            fort14 += "{}\n".format(key)
            for i in range(len(_['front_face_indexes'])):
                fort14 += "{:d} ".format(_['front_face_indexes'][i]+1)
                fort14 += "{:d} ".format(_['back_face_indexes'][i]+1)
                fort14 += "{:<.16E} ".format(_["barrier_heights"][i])
                fort14 += "{:<.16E} ".format(
                        _["subcritical_flow_coefficients"][i])
                fort14 += "{:<.16E} ".format(
                        _["supercritical_flow_coefficients"][i])
                fort14 += "\n"
        for key, _ in self.culvert_boundaries.items():
            fort14 += "{:d} ".format(len(_['indexes']))
            fort14 += "{:d} ".format(_['ibtype'])
            fort14 += "! number of nodes and ibtype for culvert_boundary_"
            fort14 += "{}\n".format(key)
            for i in range(len(_['front_face_indexes'])):
                fort14 += "{:d} ".format(_['front_face_indexes'][i]+1)
                fort14 += "{:d} ".format(_['back_face_indexes'][i]+1)
                fort14 += "{:<.16E} ".format(_["barrier_heights"][i])
                fort14 += "{:<.16E} ".format(
                        _["subcritical_flow_coefficients"][i])
                fort14 += "{:<.16E} ".format(
                        _["supercritical_flow_coefficients"][i])
                fort14 += "{:<.16E} ".format(
                        _["cross_barrier_pipe_heights"][i])
                fort14 += "{:<.16E} ".format(
                        _["friction_factors"][i])
                fort14 += "{:<.16E} ".format(
                        _["pipe_diameters"][i])
                fort14 += "\n"
        fort14 += "{}\n".format(self.SpatialReference.ExportToWkt())
        if path is not None:
            with open(path, 'w') as f:
                f.write(fort14)
        else:
            print(fort14)

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
            axes = plt.figure(figsize=figsize).add_subplot(111)
        if vmin is None:
            vmin = np.min(self.values)
        if vmax is None:
            vmax = np.max(self.values)
        cmap, norm, levels, col_val = self.__get_cmap(
            vmin, vmax, cmap, levels, colors, norm)
        axes.tricontourf(
            self.x, self.y, self.elements, self.values, levels=levels,
            cmap=cmap, norm=norm, vmin=vmin, vmax=vmax, **kwargs)
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
        cbar = plt.colorbar(mappable, cax=cax,  # extend=cmap_extend,
                            orientation='horizontal')
        if col_val != 0:
            cbar.set_ticks([vmin, vmin + col_val * (vmax-vmin), vmax])
            cbar.set_ticklabels([np.around(vmin, 2), 0.0, np.around(vmax, 2)])
        else:
            cbar.set_ticks([vmin, vmax])
            cbar.set_ticklabels([np.around(vmin, 2), np.around(vmax, 2)])
        if cbar_label is not None:
            cbar.set_label(cbar_label)
        if show is True:
            plt.show()
        return axes

    def __get_cmap(
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
                norm = _FixPointNormalize(sealevel=0.0, vmax=vmax, vmin=vmin,
                                          col_val=col_val)
        return cmap, norm, levels, col_val

    def __remove_boundary(self, key, __attribute, name):
        try:
            __attribute.pop(key)
        except KeyError:
            raise AttributeError(
                'Cannot remove {} boundary with key {}: '.format(name, key)
                + 'boundary key doesn\'t exist.')

    def __get_nodal_attribute(self, nodal_attribute):
        attribute = self.__nodal_attributes[nodal_attribute]
        attribute[nodal_attribute] = self.get_attribute(nodal_attribute)
        return attribute

    @classmethod
    def open(cls, fort14, SpatialReference=None):
        fort14 = cls.parse_fort14(fort14)
        vertices = np.vstack([fort14.pop('x'), fort14.pop('y')]).T
        cls = cls(vertices, fort14.pop('elements'), fort14.pop('values'),
                  SpatialReference, fort14.pop('node_id'),
                  fort14.pop('element_id'))
        cls.description = fort14.pop('description')
        for boundary_type in fort14.keys():
            for i, data in enumerate(fort14[boundary_type]):
                if boundary_type == 'ocean_boundaries':
                    cls.add_ocean_boundary(i, **data)
                elif boundary_type == 'land_boundaries':
                    cls.add_land_boundary(i, **data)
                elif boundary_type == 'inner_boundaries':
                    cls.add_inner_boundary(i, **data)
                elif boundary_type == 'inflow_boundaries':
                    cls.add_inflow_boundary(i, **data)
                elif boundary_type == 'outflow_boundaries':
                    cls.add_outflow_boundary(i, **data)
                elif boundary_type == 'weir_boundaries':
                    cls.add_weir_boundary(i, **data)
                elif boundary_type == 'culvert_boundaries':
                    cls.add_culvert_boundary(i, **data)
                else:
                    raise Exception('Duck-typing error.')
        return cls

    @staticmethod
    def parse_fort14(path):
        fort14 = dict()
        fort14['x'] = list()
        fort14['y'] = list()
        fort14['values'] = list()
        fort14['node_id'] = list()
        fort14['elements'] = list()
        fort14['element_id'] = list()
        fort14['ocean_boundaries'] = list()
        fort14['land_boundaries'] = list()
        fort14['inner_boundaries'] = list()
        fort14['inflow_boundaries'] = list()
        fort14['outflow_boundaries'] = list()
        fort14['weir_boundaries'] = list()
        fort14['culvert_boundaries'] = list()
        with open(path, 'r') as f:
            fort14['description'] = "{}".format(f.readline())
            NE, NP = map(int, f.readline().split())
            _NP = len([])
            while _NP < NP:
                node_id, x, y, z = f.readline().split()
                fort14['node_id'].append(int(node_id)-1)
                fort14['x'].append(float(x))
                fort14['y'].append(float(y))
                fort14['values'].append(-float(z))
                _NP += 1
            _NE = len([])
            while _NE < NE:
                line = f.readline().split()
                fort14['element_id'].append(float(line[0]))
                if int(line[1]) != 3:
                    raise NotImplementedError(
                        'Package only supports triangular meshes.')
                fort14['elements'].append([int(x)-1 for x in line[2:]])
                _NE += 1
            # Assume EOF if NOPE is empty.
            try:
                NOPE = int(f.readline().split()[0])
            except IndexError:
                return fort14
            # For now, let NOPE=-1 mean a self closing mesh
            # reassigning NOPE to 0 until further implementation is applied.
            if NOPE == -1:
                NOPE = 0
            _NOPE = len([])
            f.readline()  # Number of total open ocean nodes. Not used.
            while _NOPE < NOPE:
                fort14['ocean_boundaries'].append({'indexes': list()})
                NETA = int(f.readline().split()[0])
                _NETA = len([])
                while _NETA < NETA:
                    fort14['ocean_boundaries'][_NOPE]['indexes'].append(
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
                    fort14['land_boundaries'].append({
                                    'ibtype': IBTYPE,
                                    'indexes': list()})
                elif IBTYPE in [1, 11, 21]:
                    fort14['inner_boundaries'].append({
                                    'ibtype': IBTYPE,
                                    'indexes': list()})
                elif IBTYPE in [2, 12, 22, 102, 122]:
                    fort14['inflow_boundaries'].append({
                                    'ibtype': IBTYPE,
                                    'indexes': list()})
                elif IBTYPE in [3, 13, 23]:
                    fort14['outflow_boundaries'].append({
                                    'ibtype': IBTYPE,
                                    'indexes': list(),
                                    'barrier_heights': list(),
                                    'supercritical_flow_coefficients': list()})

                elif IBTYPE in [4, 24]:
                    fort14['weir_boundaries'].append({
                                    'ibtype': IBTYPE,
                                    'front_face_indexes': list(),
                                    'back_face_indexes': list(),
                                    'barrier_heights': list(),
                                    'subcritical_flow_coefficients': list(),
                                    'supercritical_flow_coefficients': list()})
                elif IBTYPE in [5, 25]:
                    fort14['culvert_boundaries'].append({
                                    'ibtype': IBTYPE,
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
                        fort14['land_boundaries'][-1][
                            'indexes'].append(int(line[0])-1)
                    elif IBTYPE in [1, 11, 21]:
                        fort14['inner_boundaries'][-1][
                            'indexes'].append(int(line[0])-1)
                    elif IBTYPE in [3, 13, 23]:
                        fort14['outflow_boundaries'][-1][
                            'indexes'].append(int(line[0])-1)
                        fort14['outflow_boundaries'][-1][
                            'external_barrier_heights'].append(float(line[1]))
                        fort14['outflow_boundaries'][-1][
                            'supercritical_flow_coefficients'].append(
                                float(line[2]))
                    elif IBTYPE in [2, 12, 22, 102, 122]:
                        fort14['iflowBoundaries'][-1][
                            'indexes'].append(int(line[0])-1)
                    elif IBTYPE in [4, 24]:
                        fort14['weir_boundaries'][-1][
                            'front_face_indexes'].append(int(line[0])-1)
                        fort14['weir_boundaries'][-1][
                            'back_face_indexes'].append(int(line[1])-1)
                        fort14['weir_boundaries'][-1][
                            'barrier_heights'].append(float(line[2]))
                        fort14['weir_boundaries'][-1][
                            'subcritical_flow_coefficients'].append(
                                float(line[3]))
                        fort14['weir_boundaries'][-1][
                            'supercritical_flow_coefficients'].append(
                                float(line[4]))
                    elif IBTYPE in [5, 25]:
                        fort14['culvert_boundaries'][-1][
                            'front_face_indexes'].append(int(line[0])-1)
                        fort14['culvert_boundaries'][-1][
                            'back_face_indexes'].append(int(line[1])-1)
                        fort14['culvert_boundaries'][-1][
                            'barrier_heights'].append(float(line[2]))
                        fort14['culvert_boundaries'][-1][
                            'subcritical_flow_coefficients'].append(
                                float(line[3]))
                        fort14['culvert_boundaries'][-1][
                            'supercritical_flow_coefficients'].append(
                                float(line[4]))
                        fort14['culvert_boundaries'][-1][
                            'friction_factors'].append(float(line[5]))
                        fort14['culvert_boundaries'][-1][
                            'pipe_diameters'].append(float(line[6]))
                    else:
                        Exception("Duck-typing error. "
                                  + "This exception should be unreachable.")
                    _NVELL += 1
                _NBOU += 1
        return fort14

    @staticmethod
    def parse_fort13(path):
        fort13 = dict()
        with open(path, 'r') as f:
            f.readline()
            f.readline()
            NAttr = int(f.readline().split()[0])
            i = 0
            while i < NAttr:
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
            for i in range(NAttr):
                attribute_name = f.readline().strip()
                numOfNodes = int(f.readline())
                indexes = list()
                values = list()
                j = 0
                while j < numOfNodes:
                    line = f.readline().split()
                    indexes.append(int(line.pop(0))-1)
                    line = [float(x) for x in line]
                    values.append(line)
                    j += 1
                fort13[attribute_name]['indexes'] = indexes
                fort13[attribute_name]['values'] = values
        return fort13

    @property
    def description(self):
        if not hasattr(self, "__description"):
            self.__description = "AdcircPy mesh created on " \
                                 + "{}".format(datetime.now())
        return self.__description

    @property
    def node_id(self):
        if self.__node_id is not None:
            return self.__node_id
        else:
            node_id = np.arange(self.vertices.shape[0])
            self.__node_id = np.asarray(node_id, dtype=np.int64)
            return self.__node_id

    @property
    def element_id(self):
        if self.__element_id is not None:
            return self.__element_id
        else:
            element_id = np.arange(self.elements.shape[0])
            self.__element_id = np.asarray(element_id, dtype=np.int64)
            return self.__element_id

    @property
    def ocean_boundary(self):
        ocean_boundary = list()
        for boundary in self.ocean_boundaries.values():
            for idx in boundary:
                ocean_boundary.append(idx)
        return ocean_boundary

    @property
    def land_boundary(self):
        land_boundary = list()
        for boundary in self.land_boundaries.values():
            for idx in boundary['indexes']:
                land_boundary.append(idx)
        return land_boundary

    @property
    def inner_boundary(self):
        inner_boundary = list()
        for key, boundary in self.inner_boundaries.items():
            for idx in boundary:
                inner_boundary.append(idx)
        return inner_boundary

    @property
    def inflow_boundary(self):
        inflow_boundary = list()
        for key, boundary in self.inflow_boundaries.items():
            for idx in boundary:
                inflow_boundary.append(idx)
        return inflow_boundary

    @property
    def outflow_boundary(self):
        outflow_boundary = list()
        for key, boundary in self.outflow_boundaries.items():
            for idx in boundary:
                outflow_boundary.append(idx)
        return outflow_boundary

    @property
    def weir_boundary(self):
        weir_boundary = list()
        for key, boundary in self.weir_boundaries.items():
            for idx in boundary:
                weir_boundary.append(idx)
        return weir_boundary

    @property
    def culvert_boundary(self):
        culvert_boundary = list()
        for key, boundary in self.culvert_boundaries.items():
            for idx in boundary:
                culvert_boundary.append(idx)
        return culvert_boundary

    @property
    def ocean_boundaries(self):
        ocean_boundaries = dict()
        __ocean_boundaries = deepcopy(self.__ocean_boundaries)
        for key, _ in __ocean_boundaries.items():
            ocean_boundaries[key] = list(np.where(_.pop('__bool'))[0])
        return ocean_boundaries

    @property
    def land_boundaries(self):
        land_boundaries = dict()
        __land_boundaries = deepcopy(self.__land_boundaries)
        for key, _ in __land_boundaries.items():
            land_boundaries[key] = {
                'indexes': list(np.where(_.pop('__bool'))[0]), **_}
        return land_boundaries

    @property
    def inner_boundaries(self):
        inner_boundaries = dict()
        __inner_boundaries = deepcopy(self.__inner_boundaries)
        for key, _ in __inner_boundaries.items():
            inner_boundaries[key] = {
                'indexes': list(np.where(_.pop('__bool'))[0]), **_}
        return inner_boundaries

    @property
    def inflow_boundaries(self):
        inflow_boundaries = dict()
        __inflow_boundaries = deepcopy(self.__inflow_boundaries)
        for key, _ in __inflow_boundaries.items():
            inflow_boundaries[key] = {
                'indexes': list(np.where(_.pop('__bool'))[0]), **_}
        return inflow_boundaries

    @property
    def outflow_boundaries(self):
        outlow_boundaries = dict()
        __outflow_boundaries = deepcopy(self.__outflow_boundaries)
        for key, _ in __outflow_boundaries.items():
            outlow_boundaries[key] = {
                'indexes': list(np.where(_.pop('__bool'))[0]), **_}
        return outlow_boundaries

    @property
    def weir_boundaries(self):
        weir_boundaries = dict()
        __weir_boundaries = deepcopy(self.__weir_boundaries)
        for key, _ in __weir_boundaries.items():
            front_face = list(np.where(_.pop('__bool_front_face'))[0])
            back_face = list(np.where(_.pop('__bool_back_face'))[0])
            weir_boundaries[key] = {'front_face_indexes': front_face,
                                    'back_face_indexes': back_face, **_}
        return weir_boundaries

    @property
    def culvert_boundaries(self):
        culvert_boundaries = dict()
        __culvert_boundaries = deepcopy(self.__culvert_boundaries)
        for key, _ in __culvert_boundaries.items():
            front_face = list(np.where(_.pop('__bool_front_face'))[0])
            back_face = list(np.where(_.pop('__bool_back_face'))[0])
            culvert_boundaries[key] = {'front_face_indexes': front_face,
                                       'back_face_indexes': back_face, **_}
        return culvert_boundaries

    @property
    def primitive_weighting_in_continuity_equation(self):
        return self.__get_nodal_attribute(
            "primitive_weighting_in_continuity_equation")

    @property
    def surface_submergence_state(self):
        return self.__get_nodal_attribute("surface_submergence_state")

    @property
    def quadratic_friction_coefficient_at_sea_floor(self):
        return self.__get_nodal_attribute(
            "quadratic_friction_coefficient_at_sea_floor")

    @property
    def surface_directional_effective_roughness_length(self):
        return self.__get_nodal_attribute(
            "surface_directional_effective_roughness_length")

    @property
    def surface_canopy_coefficient(self):
        return self.__get_nodal_attribute("surface_canopy_coefficient")

    @property
    def bridge_pilings_friction_parameters(self):
        return self.__get_nodal_attribute("bridge_pilings_friction_parameters")

    @property
    def mannings_n_at_sea_floor(self):
        return self.__get_nodal_attribute("mannings_n_at_sea_floor")

    @property
    def chezy_friction_coefficient_at_sea_floor(self):
        return self.__get_nodal_attribute(
            "chezy_friction_coefficient_at_sea_floor")

    @property
    def sea_surface_height_above_geoid(self):
        return self.__get_nodal_attribute("sea_surface_height_above_geoid")

    @property
    def bottom_roughness_length(self):
        return self.__get_nodal_attribute("bottom_roughness_length")

    @property
    def wave_refraction_in_swan(self):
        return self.__get_nodal_attribute("wave_refraction_in_swan")

    @property
    def average_horizontal_eddy_viscosity_in_sea_water_wrt_depth(self):
        return self.__get_nodal_attribute(
            "average_horizontal_eddy_viscosity_in_sea_water_wrt_depth")

    @property
    def elemental_slope_limiter(self):
        return self.__get_nodal_attribute("elemental_slope_limiter")

    @property
    def advection_state(self):
        return self.__get_nodal_attribute("advection_state")

    @property
    def initial_river_elevation(self):
        return self.__get_nodal_attribute("initial_river_elevation")

    @description.setter
    def description(self, description):
        self.__decription = str(description)

    @node_id.setter
    def node_id(self, node_id):
        if node_id is not None:
            node_id = np.asarray(node_id)
            assert node_id.shape[0] == self.vertices.shape[0]
        self.__node_id = node_id

    @element_id.setter
    def element_id(self, element_id):
        if element_id is not None:
            element_id = np.asarray(element_id)
            assert element_id.shape[0] == self.elements.shape[0]
        self.__element_id = element_id
