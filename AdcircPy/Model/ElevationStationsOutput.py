from datetime import timedelta
from collections import OrderedDict
from AdcircPy.Model._StationsOutput import _StationsOutput


class ElevationStationsOutput(_StationsOutput):

    __storage = OrderedDict()

    def __init__(self, sampling_frequency=timedelta(0), netcdf=True,
                 spinup=False, harmonic_analysis=False):
        super(ElevationStationsOutput, self).__init__(
            sampling_frequency, netcdf, spinup, harmonic_analysis)

    def __getitem__(self, key):
        return self._storage[key]

    def __iter__(self):
        return iter(self._storage)

    def __len__(self):
        return len(self._storage.keys())

    def add_stations_from_fort15(self, path, _hint='NOUTE'):
        super(ElevationStationsOutput, self).add_stations_from_fort15(path,
                                                                      _hint)

    @classmethod
    def from_fort15(cls, path, sampling_frequency, netcdf=True,
                    spinup=False, harmonic_analysis=False, _hint='NOUTE'):
        return cls._from_fort15(path, _hint, sampling_frequency, netcdf,
                                spinup, harmonic_analysis)

    @property
    def _storage(self):
        return self.__storage
