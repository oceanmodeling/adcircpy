import logging
import pathlib

import numpy as np

logger = logging.getLogger(__name__)


class NodalAttributes:
    def __init__(self, fort14):
        self._fort14 = fort14
        self._attributes = {}

    def __iter__(self):
        for name, data in self._attributes.items():
            yield name, data

    def __str__(self):
        fort13 = [
            f'{self._fort14.description} nodal attributes',
            f'{len(self._fort14.nodes.id)}',
            f'{len(self.get_attribute_names())}',
        ]

        for name in self.get_attribute_names():
            attribute = self.get_attribute(name)
            fort13.extend(
                [
                    f'{name}',
                    f'{attribute["units"]}',
                    f'{len(attribute["defaults"])}',
                    ' '.join(f'{n:<.16E}' for n in attribute['defaults']),
                ]
            )
        for name in self.get_attribute_names():
            attribute = self.get_attribute(name)
            fort13.extend(
                [f'{name}', f'{len(attribute["non_default_indexes"])}',]
            )
            for i, values in enumerate(
                attribute['values'][attribute['non_default_indexes'], :]
            ):
                node_id = self._fort14.nodes.get_id_by_index(
                    attribute['non_default_indexes'][i]
                )
                line = [f'{node_id}']
                for value in values:
                    line.append(f'{value}')
                fort13.append(' '.join(line))
        return '\n'.join(fort13)

    def add_attribute(self, name: str, units: str = None):
        if name in self._attributes:
            raise AttributeError(
                f'Cannot add nodal attribute with name ' f'{name}: attribute already exists.'
            )
        self._attributes.setdefault(name, {}).update(
            {
                'units': 'unitless' if units is None else str(units),
                'values': None,
                'coldstart': False,
                'hotstart': False,
            }
        )

    def get_coldstart_attributes(self):
        coldstart_attributes = dict()
        for attribute in self.get_attribute_names():
            attr = self.get_attribute(attribute)
            if attr['coldstart']:
                coldstart_attributes[attribute] = attr
        return coldstart_attributes

    def get_hotstart_attributes(self):
        hotstart_attributes = dict()
        for attribute in self.get_attribute_names():
            attr = self.get_attribute(attribute)
            if attr['hotstart']:
                hotstart_attributes[attribute] = attr
        return hotstart_attributes

    def set_attribute_coldstart_state(self, attribute, state):
        assert isinstance(state, bool)
        self.get_attribute(attribute)['coldstart'] = state

    def set_attribute_hotstart_state(self, attribute, state):
        assert isinstance(state, bool)
        self.get_attribute(attribute)['hotstart'] = state

    def set_attribute_state(self, attribute, coldstart, hotstart):
        self.set_attribute_coldstart_state(attribute, coldstart)
        self.set_attribute_hotstart_state(attribute, hotstart)

    def get_attribute_names(self):
        return self._attributes.keys()

    def get_attribute(self, name):
        if name not in self.get_attribute_names():
            msg = f'Nodal attribute with name {name} has not been loaded.'
            raise AttributeError(msg)
        if self._attributes[name].get('defaults') is None:
            # TODO: the 'values' array can be generated more succintly.
            def mode_rows(a):
                a = np.ascontiguousarray(a)
                void_dt = np.dtype((np.void, a.dtype.itemsize * np.prod(a.shape[1:])))
                _, ids, count = np.unique(
                    a.view(void_dt).ravel(), return_index=True, return_counts=True
                )
                largest_count_id = ids[count.argmax()]
                most_frequent_row = a[largest_count_id]
                return most_frequent_row

            _attr = self._attributes[name]
            if _attr['values'].ndim == 1:  # rewrite as column major
                _attr['values'] = _attr['values'].reshape((_attr['values'].shape[0], 1))
            _attr['defaults'] = mode_rows(_attr['values'])
            _attr['non_default_indexes'] = np.where(
                (_attr['values'] != _attr['defaults']).all(axis=1)
            )[0]
            self._attributes[name] = _attr
        return self._attributes[name]

    def has_attribute(self, attribute_name, runtype=None):
        """
        if runtype is None, returns True if attribute_name is in any.
        """
        if attribute_name not in self.get_attribute_names():
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

    def set_attribute(self, name, values, coldstart: bool = False, hotstart: bool = False):
        if name not in self.get_attribute_names():
            raise AttributeError(
                f'Cannot set nodal attribute with name '
                f'{name}: attribute has not been '
                f'added yet.'
            )
        assert isinstance(coldstart, bool)
        assert isinstance(hotstart, bool)
        assert values.flatten().shape[0] % self._fort14.values.shape[0] == 0
        self._attributes[name].update(
            {'values': values, 'coldstart': coldstart, 'hotstart': hotstart}
        )

    def add_patch(self, name, patch, value):
        if name not in self.get_attribute_names():
            raise AttributeError(
                f'Cannot add patch to nodal attribute with name {name}: '
                'attribute has not been added yet.'
            )
        idxs = [
            row.Index for row in self._fort14.nodes.gdf.geometry.within(patch).itertuples()
        ]
        self._attributes[name]['values'][idxs] = value

    def import_fort13(self, fort13):
        fort13 = parse_fort13(fort13)
        if fort13.pop('NumOfNodes') != len(self._fort14.nodes.id):
            raise Exception('fort.13 file does not match the mesh.')
        self._AGRID = fort13.pop('AGRID')
        for attribute, data in fort13.items():
            values = np.asarray(data['values'])
            if values.ndim == 1:
                values = values.reshape((values.shape[0], 1))
            full_values = np.full(
                (self._fort14.values.size, np.asarray(data['defaults']).flatten().size),
                np.nan,
            )
            for i, idx in enumerate(data['indexes']):
                for j, value in enumerate(values[i, :].tolist()):
                    full_values[idx, j] = value
            idxs = np.where(np.isnan(full_values).all(axis=1))[0]
            for idx in idxs:
                for i, value in enumerate(data['defaults']):
                    full_values[idx, i] = value
            # converts from column major to row major, leave it column major.
            # if full_values.shape[1] == 1:
            #     full_values = full_values.flatten()
            self.add_attribute(attribute, data['units'])
            self.set_attribute(attribute, full_values)

    def write(self, path, overwrite=False):
        if path is not None:
            path = pathlib.Path(path)
            if path.is_file() and not overwrite:
                msg = 'File exists, pass overwrite=True to allow overwrite.'
                raise Exception(msg)
            else:
                with open(path, 'w') as f:
                    f.write(str(self))
        else:
            print(str(self))


def parse_fort13(path):
    fort13 = {}
    with open(path, 'r') as f:
        fort13['AGRID'] = f.readline().strip()
        fort13['NumOfNodes'] = NP = int(f.readline().split()[0])
        NAttr = int(f.readline().split()[0])
        i = 0
        while i < NAttr:
            attribute_name = f.readline().strip()
            units = f.readline().strip()
            if units == '1':
                units = 'unitless'
            f.readline()
            defaults = [float(x) for x in f.readline().split()]
            fort13[attribute_name] = {
                'units': units,
                'defaults': defaults,
                'indexes': [],
            }
            i += 1
        for i in range(NAttr):
            attribute_name = f.readline().strip()
            numOfNodes = int(f.readline())
            values = np.zeros((NP, len(fort13[attribute_name]['defaults'])))
            values[:] = np.nan
            j = 0
            while j < numOfNodes:
                str = f.readline().split()
                index = int(str[0]) - 1
                fort13[attribute_name]['indexes'].append(index)
                node_values = [float(x) for x in str[1:]]
                values[index, :] = node_values
                j += 1
            values[np.where(np.isnan(values[:, 0])), :] = fort13[attribute_name]['defaults']
            fort13[attribute_name]['values'] = values
        return fort13
