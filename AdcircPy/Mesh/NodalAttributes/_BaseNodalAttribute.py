# global imports
import numpy as np
from scipy import stats


class _BaseNodalAttribute(object):

    def __init__(self, UnstructuredMesh, values, default_values=None,
                 units=None):
        values = np.asarray(values)
        assert UnstructuredMesh.xy.shape[0] == values.shape[0]
        self._values = values
        self._default_values = default_values
        self._units = units
        self.__set_indexes()

    def __set_indexes(self):
        indexes = np.where((self.values != self.default_values).all(axis=1))[0]
        self.__indexes = indexes

    @classmethod
    def _from_fort13(cls, UnstructuredMesh, indexes, values, default_values,
                     units):
        default_values = np.asarray(default_values).flatten()
        full_values = np.full((UnstructuredMesh.xy.shape[0],
                               default_values.shape[0]), np.nan)
        full_values[indexes, :] = values
        nan_indexes = np.where(np.isnan(full_values))
        for j, default in enumerate(default_values):
            full_values[nan_indexes, j] = default
        return cls(UnstructuredMesh, values=full_values,
                   default_values=default_values, units=units)

    @property
    def values(self):
        return self._values

    @property
    def default_values(self):
        return self._default_values

    @property
    def units(self):
        return self._units

    @property
    def indexes(self):
        return self.__indexes

    @property
    def _indexes(self):
        return self.__indexes

    @property
    def _values(self):
        return self.__values

    @property
    def _default_values(self):
        return self.__default_values

    @property
    def _units(self):
        return self.__units

    @_indexes.setter
    def _indexes(self, indexes):
        if indexes is not None:
            raise NotImplementedError
        self.__indexes = indexes

    @_values.setter
    def _values(self, values):
        if values is None:
            raise NotImplementedError
        values = np.asarray(values)
        if values.shape[0] == values.size:
            values = values.reshape((values.size, 1))
        self.__values = values

    @_default_values.setter
    def _default_values(self, default_values):
        if default_values is None:
            raise NotImplementedError('Must pass vector of default values for '
                                      + 'now.')
            default_values = stats.mode(self.values)

        self.__default_values = np.asarray(default_values).flatten()

    @_units.setter
    def _units(self, units):
        self.__units = units
