from datetime import timedelta
from AdcircPy.Model._StationsOutput import _StationsOutput


class VelocityStationsOutput(_StationsOutput):

    def __init__(self, sampling_frequency=timedelta(0),
                 netcdf=True, spinup=False, harmonic_analysis=False):
        self.__storage = dict()
        super(VelocityStationsOutput, self).__init__(
            sampling_frequency, netcdf, spinup, harmonic_analysis)

    def add_stations_from_fort15(self, path, _hint='NOUTV'):
        super(VelocityStationsOutput, self).add_stations_from_fort15(path,
                                                                     _hint)

    @classmethod
    def from_fort15(cls, path, sampling_frequency=None, netcdf=True,
                    spinup=False, harmonic_analysis=False, _hint='NOUTV'):
        return super(VelocityStationsOutput, cls).__from_fort15(
            path, _hint, sampling_frequency, netcdf, spinup,
            harmonic_analysis)

    @property
    def _storage(self):
        return self.__storage
