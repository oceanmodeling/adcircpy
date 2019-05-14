from datetime import timedelta
from collections import OrderedDict
from AdcircPy.Model._StationsOutput import _StationsOutput


class MeteorologicalStationsOutput(_StationsOutput):

    __storage = OrderedDict()

    def __init__(self, sampling_frequency=timedelta(0), netcdf=True):
        super(MeteorologicalStationsOutput, self).__init__(
                                    sampling_frequency, netcdf, False, False)

    def __getitem__(self, key):
        return self._storage[key]

    def __iter__(self):
        return iter(self._storage)

    def __len__(self):
        return len(self._storage.keys())

    def add_stations_from_fort15(self, path, _hint='NOUTM'):
        super(MeteorologicalStationsOutput, self).add_stations_from_fort15(
                                                                path, _hint)

    @classmethod
    def from_fort15(cls, path, sampling_frequency=None, netcdf=True,
                    spinup=False, harmonic_analysis=False, _hint='NOUTV'):
        return super(MeteorologicalStationsOutput, cls).__from_fort15(
            path, _hint, sampling_frequency, netcdf, spinup,
            harmonic_analysis)

    @property
    def _storage(self):
        return self.__storage
