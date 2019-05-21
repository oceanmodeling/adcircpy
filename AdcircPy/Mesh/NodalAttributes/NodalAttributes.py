# global imports
from collections.abc import Mapping
# global imports
from pathlib import Path

# local imports
from AdcircPy.Mesh.NodalAttributes import _PrimitiveWeighting
from AdcircPy.Mesh.NodalAttributes import _BaseNodalAttribute


class NodalAttributes(Mapping):

    __allowed_attributes = [
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
        'initial_river_elevation']
    __storage = dict()
    __coldstart_attributes = set()
    __hotstart_attributes = set()
    __all_attributes = set()

    def __init__(self, UnstructuredMesh, coldstart_attributes=[],
                 hotstart_attributes=[], use_all=False, **nodal_attributes):
        # assert isinstance(UnstructuredMesh, _UnstructuredMesh)
        self.AGRID = UnstructuredMesh.description
        self.NumOfNodes = UnstructuredMesh.xy.shape[0]
        for name, values in nodal_attributes.items():
            self.add_attribute(UnstructuredMesh, name, **values)
            if use_all is True:
                self.set_attribute_state(name, coldstart=True, hotstart=True)
            else:
                if name in coldstart_attributes:
                    coldstart = True
                else:
                    coldstart = False
                if name in hotstart_attributes:
                    hotstart = True
                else:
                    hotstart = False
                self.set_attribute_state(name, coldstart=coldstart,
                                         hotstart=hotstart)

    def __getitem__(self, key):
        return self._storage[key]

    def __iter__(self):
        return iter(self._storage.items())

    def __len__(self):
        return len(self._storage.keys())

    def __repr__(self):
        print('{}'.format(self._storage))

    def set_attribute_state(self, attribute, coldstart, hotstart):
        if attribute not in self._storage.keys():
            raise RuntimeError('set_attribute_state(): {} '.format(attribute)
                               + 'is not a loaded attribute.')
        if coldstart is True:
            if attribute not in self.coldstart_attributes:
                self.coldstart_attributes.add(attribute)
            else:
                self.coldstart_attributes.remove(attribute)
        else:
            if attribute in self.coldstart_attributes:
                self.coldstart_attributes.remove(attribute)
            else:
                self.coldstart_attributes.add(attribute)

        if hotstart is True:
            if attribute not in self.hotstart_attributes:
                self.hotstart_attributes.add(attribute)
            else:
                self.hotstart_attributes.remove(attribute)
        else:
            if attribute in self.hotstart_attributes:
                self.hotstart_attributes.remove(attribute)
            else:
                self.hotstart_attributes.add(attribute)

    def set_TAU0(self, UnstructuredMesh, TAU0, **kwargs):
        attr_name = 'primitive_weighting_in_continuity_equation'
        if attr_name in self._storage.keys():
            if TAU0 == 3.:
                self._storage[attr_name] = _PrimitiveWeighting.from_TAU0(
                    UnstructuredMesh, TAU0, **kwargs)
        else:
            if TAU0 is not None:
                self._storage[attr_name] = _PrimitiveWeighting.from_TAU0(
                    UnstructuredMesh, TAU0, **kwargs)

    def auto_generate_TAU0(self, **kwargs):
        self._storage['primitive_weighting_in_continuity_equation'] \
            = _PrimitiveWeighting(self, 3, **kwargs)

    def dump(self, output_dir=None, filename='fort.13'):
        if output_dir is not None:
            output_path = Path(str(output_dir) + '/' + filename)
            if len(self._storage.keys()) > 0:
                fort13 = self.get_fort13()
                with open(str(output_path), 'w') as f:
                    f.write(fort13)
        else:
            if len(self._storage.keys()) > 0:
                print(self.get_fort13())

    def add_attribute(self, UnstructuredMesh, attribute_name,
                      **attribute_values):
        if attribute_name == 'primitive_weighting_in_continuity_equation':
            self._storage[attribute_name] = _PrimitiveWeighting.from_TAU0(
                                    UnstructuredMesh, -3, **attribute_values)
        elif attribute_name in self.allowed_attributes:
            self._storage[attribute_name] = _BaseNodalAttribute._from_fort13(
                                        UnstructuredMesh, **attribute_values)
        else:
            raise AttributeError('Attribute '
                                 + '{} '.format(attribute_name)
                                 + 'is not a known ADCIRC attribute')

    def get_fort13(self):
        f = ''
        f += "{}\n".format(self.AGRID)
        f += "{}\n".format(self.NumOfNodes)
        NAttr = len(self._storage.keys())
        f += "{}\n".format(NAttr)
        for name, attribute in self:
            f += "{}\n".format(name)
            f += "{}\n".format(attribute.units)
            f += "{}\n".format(len(attribute.default_values))
            for value in attribute.default_values:
                f += "{} ".format(value)
            f += "\n"
        for name, attribute in self:
            f += "{}\n".format(name)
            f += "{}\n".format(len(attribute.indexes))
            for i, idx in enumerate(attribute.indexes):
                f += "{:<10d} ".format(idx+1)
                for value in attribute.values[idx, :]:
                    f += "{} ".format(value)
                f += "\n"
        return f

    @classmethod
    def from_fort13(cls, UnstructuredMesh, path, coldstart_attributes=[],
                    hotstart_attributes=[], use_all=False):
        coldstart_attributes = list(coldstart_attributes)
        hotstart_attributes = list(hotstart_attributes)
        return cls(UnstructuredMesh, coldstart_attributes, hotstart_attributes,
                   use_all, **cls.parse_fort13(path))

    @staticmethod
    def parse_fort13(path):
        data = dict()
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
                data[attribute_name] = dict()
                data[attribute_name]['units'] = units
                data[attribute_name]['default_values'] = defaults
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
                data[attribute_name]['indexes'] = indexes
                data[attribute_name]['values'] = values
        return data

    @property
    def coldstart_attributes(self):
        """ """
        return self.__coldstart_attributes

    @property
    def hotstart_attributes(self):
        """ """
        return self.__hotstart_attributes

    @property
    def allowed_attributes(self):
        """ """
        return self.__allowed_attributes

    @property
    def coldstart(self):
        """ """
        return self.__coldstart_attributes

    @property
    def hotstart(self):
        """ """
        return self.__hotstart_attributes

    @property
    def TAU0(self):
        """ """
        if 'primitive_weighting_in_continuity_equation' in self:
            return self._storage[
                'primitive_weighting_in_continuity_equation'].TAU0
        else:
            return None

    @property
    def _storage(self):
        return self.__storage
