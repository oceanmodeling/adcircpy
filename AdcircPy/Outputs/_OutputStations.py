from collections.abc import Mapping
from AdcircPy.Outputs._StationTimeseries import _StationTimeseries


class _OutputStations(Mapping):

    def __init__(self, **stations):
        self.__set_storage(stations)

    def __getitem__(self, key):
        if isinstance(key, int):
            for i, _key in enumerate(self._storage.keys()):
                if i == key:
                    return self._storage[_key]
        else:
            return self._storage[key]

    def __iter__(self):
        return iter(self._storage)

    def __len__(self):
        return len(self._storage.keys())

    def __set_storage(self, stations):
        self.__storage = dict()
        for station_id, data in stations.items():
            self._storage[station_id] = _StationTimeseries(**data)

    @property
    def _storage(self):
        return self.__storage
