# type: ignore[attr-defined]
from collections import defaultdict
from itertools import permutations

from haversine import haversine, Unit
import numpy as np

from adcircpy.forcing.bctypes import BoundaryCondition
from adcircpy.forcing.tides import Tides
from adcircpy.forcing.waves import WaveForcing
from adcircpy.forcing.winds.base import WindForcing
from adcircpy.mesh.fort13 import NodalAttributes
from adcircpy.mesh.fort14 import Fort14


class ModelForcings:
    def __init__(self, fort14):
        self.wind = None
        self.wave = None
        self.tides = None

    def add(self, forcing):
        if isinstance(forcing, BoundaryCondition):
            if isinstance(forcing, Tides):
                self.tides = forcing
            else:
                raise NotImplementedError(
                    f'Unhandled boundary condition of type {type(forcing)}.'
                )

        elif isinstance(forcing, WindForcing):
            self.wind = forcing

        elif isinstance(forcing, WaveForcing):
            self.wave = forcing

        else:
            msg = f'Unrecognized forcing type {forcing}.'
            raise Exception(msg)

    def __eq__(self, other: 'ModelForcings') -> bool:
        return self.__class__ == other.__class__ and self.__dict__ == other.__dict__


class NodalAttributeDescriptor:
    def __init__(self, name):
        self.name = name

    def __set__(self, obj, val):
        if not obj.nodal_attributes.has_attribute(self.name):
            obj.nodal_attributes.add_attribute(self.name)
        obj.nodal_attributes.set_attribute(self.name, val, True, True)

    def __get__(self, obj, val):
        return self.nodal_attributes.get_attribute(self.name)


class AdcircMeshMeta(type):
    adcirc_nodal_attributes = [
        'primitive_weighting_in_continuity_equation',
        'surface_submergence_state',
        'quadratic_friction_coefficient_at_sea_floor',
        'surface_directional_effective_roughness_length',
        'surface_canopy_coefficient',
        'bridge_pilings_friction_parameters',
        'mannings_n_at_sea_floor',
        'chezy_friction_coefficient_at_sea_floor',
        'sea_surface_height_above_geoid',
        'bottom_roughness_length',
        'wave_refraction_in_swan',
        'average_horizontal_eddy_viscosity_in_sea_water_wrt_depth',
        'elemental_slope_limiter',
        'advection_state',
        'initial_river_elevation',
    ]

    def __new__(meta, name, bases, attrs):
        for attribute in meta.adcirc_nodal_attributes:
            attrs[attribute] = NodalAttributeDescriptor(attribute)
        return type(name, (Fort14,), attrs)


class AdcircMesh(metaclass=AdcircMeshMeta):
    """
    Class used to configure ADCIRC model runs.
    """

    def __init__(self, *args, fort13=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.forcings = ModelForcings(self)
        self.nodal_attributes = NodalAttributes(self)
        if fort13 is not None:
            self.import_nodal_attributes(fort13)

    def add_forcing(self, forcing):
        self.forcings.add(forcing)

    def add_nodal_attribute(self, name: str, units: str):
        self.nodal_attributes.add_attribute(name, units)

    def set_nodal_attribute(
        self, name, values, coldstart: bool = False, hotstart: bool = False
    ):
        self.nodal_attributes.set_attribute(name, values, coldstart, hotstart)

    def get_coldstart_nodal_attributes(self):
        return self.nodal_attributes.get_coldstart_attributes()

    def get_hotstart_nodal_attributes(self):
        return self.nodal_attributes.get_hotstart_attributes()

    def set_nodal_attribute_coldstart_state(self, name, state):
        self.nodal_attributes.set_attribute_coldstart_state(name, state)

    def set_nodal_attribute_hotstart_state(self, name, state):
        self.nodal_attributes.set_attribute_hotstart_state(name, state)

    def set_nodal_attribute_state(self, name, coldstart, hotstart):
        self.nodal_attributes.set_attribute_state(name, coldstart, hotstart)

    def get_nodal_attribute_names(self):
        return self.nodal_attributes.get_attribute_names()

    def get_nodal_attribute(self, name):
        return self.nodal_attributes.get_attribute(name)

    def add_nodal_attribute_patch(self, name, patch, value):
        self.nodal_attributes.add_patch(name, patch, value)

    def has_nodal_attribute(self, name, runtype=None):
        return self.nodal_attributes.has_attribute(name, runtype)

    def import_nodal_attributes(self, fort13, enable: bool = False):
        self.nodal_attributes.import_fort13(fort13)
        if bool(enable) is True:
            for attribute in self.get_nodal_attribute_names():
                self.set_nodal_attribute_state(attribute, True, True)

    def generate_constant_mannings_n(self, value: float):
        self.mannings_n_at_sea_floor = self.coords.shape[0] * [value]

    def generate_linear_mannings_n(
        self,
        min_value: float = 0.02,
        max_value: float = 0.05,
        min_depth: float = None,
        max_depth: float = None,
    ):

        # Inspired by https://github.com/schism-dev/schism/blob/master/src/Utility/Pre-Processing/NWM/Manning/write_manning.py

        min_depth = np.min(self.values) if min_depth is None else float(min_depth)
        max_depth = np.max(self.values) if max_depth is None else float(max_depth)

        values = min_value + (self.values - min_depth) * (max_value - min_value) / (
            max_depth - min_depth
        )

        if min_value is not None:
            values[values < min_value] = min_value

        if max_value is not None:
            values[values > max_value] = max_value

        self.mannings_n_at_sea_floor = values

    def generate_tau0(
        self,
        default_value=0.03,
        threshold_distance=1750.0,
        shallow_tau0=0.02,
        deep_tau0=0.005,
        threshold_depth=-10.0,
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
        msg = 'Cannot compute TAU0 with nan depth values.'
        assert not np.any(np.isnan(self.values)), msg
        msg = 'Cannot compute TAU0 with no coordinate reference system set.'
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
                if self.values.iloc[k, :].values[0] >= threshold_depth:
                    values[k] = shallow_tau0
                else:
                    values[k] = deep_tau0
        self.primitive_weighting_in_continuity_equation = values

    def critical_timestep(self, cfl, maxvel=5.0, g=9.8):
        """
        http://swash.sourceforge.net/online_doc/swashuse/node47.html
        """
        dxdy = len(self.values) * [None]
        for k, v in self.node_distances_in_meters.items():
            _dxdy = []
            for idx in v:
                _dxdy.append(self.node_distances_in_meters[k][idx])
            dxdy[k] = np.min(_dxdy)
        return cfl * np.min(dxdy) / np.abs(maxvel)

    @property
    def node_distances_in_meters(self):
        if not hasattr(self, '_node_distances_in_meters'):
            points = self.get_xy('EPSG:4326')
            self._node_distances_in_meters = {}
            for k, v in self.node_neighbors.items():
                x0, y0 = points.iloc[k].values
                self._node_distances_in_meters[k] = {}
                for idx in v:
                    x1, y1 = points.iloc[idx].values
                    self._node_distances_in_meters[k][idx] = haversine(
                        (x0, y0), (x1, y1), unit=Unit.METERS
                    )
        return self._node_distances_in_meters

    @property
    def node_neighbors(self):
        if not hasattr(self, '_node_neighbors'):
            self._node_neighbors = defaultdict(set)
            for simplex in self.triangulation.triangles:
                for i, j in permutations(simplex, 2):
                    self._node_neighbors[i].add(j)
        return self._node_neighbors

    def __copy__(self) -> bool:
        instance = super().__copy__()
        instance.forcings = self.forcings
        instance.nodal_attributes = self.nodal_attributes
        return instance

    def __eq__(self, other: 'AdcircMesh') -> bool:
        return (
            super().__eq__(other)
            and self.forcings == other.forcings
            and self.nodal_attributes == other.nodal_attributes
        )
