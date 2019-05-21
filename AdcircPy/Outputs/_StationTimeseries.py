import numpy as np


class _StationTimeseries(object):

    def __init__(self, x, y, values, time, name=None):
        self._x = x
        self._y = y
        self._values = values
        self._time = time
        assert self.values.size == self.time.size
        self._name = name

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def values(self):
        return self._values

    @property
    def time(self):
        return self._time

    @property
    def name(self):
        return self._name

    @property
    def _station_id(self):
        return self.__station_id

    @property
    def _x(self):
        return self.__x

    @property
    def _y(self):
        return self.__y

    @property
    def _values(self):
        return self.__values

    @property
    def _time(self):
        return self.__time

    @property
    def _name(self):
        return self.__name

    @_x.setter
    def _x(self, x):
        self.__x = float(x)

    @_y.setter
    def _y(self, y):
        self.__y = float(y)

    @_values.setter
    def _values(self, values):
        self.__values = np.asarray(values).flatten()

    @_time.setter
    def _time(self, time):
        self.__time = np.asarray(time).flatten()

    @_name.setter
    def _name(self, name):
        self.__name = str(name).strip('\n')
