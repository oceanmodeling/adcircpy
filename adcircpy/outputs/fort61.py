from datetime import datetime, timedelta
import pathlib
import uuid

from netCDF4 import Dataset

"""
This class could probably be improved by following the suggestions found
https://stackoverflow.com/questions/4014621/a-python-class-that-acts-like-dict
"""


class ElevationStations:

    def __init__(self, path):
        self._path = path

    def __iter__(self):
        for key, item in self.stations.items():
            yield key, item

    def __repr__(self):
        print(self.stations)

    def _certify_netcdf_stations_file(self, nc):
        msg = f'Input file {self.path} is not an ADCIRC stations output file '
        msg += '(fort.61.nc).'
        assert 'station_name' in nc.variables and 'zeta' in nc.variables, msg

    def _init_netcdf_stations(self):
        stations = dict()
        for idx, name in enumerate(self.nc['station_name']):
            name = ''.join(
                [s.decode('UTF-8') for s in name]).strip(' ')
            if len(name) == 0 or name in stations.keys():
                name = uuid.uuid4().hex[:8]
            stations[name] = dict()
            stations[name]['x'] = float(self.nc['x'][idx])
            stations[name]['y'] = float(self.nc['y'][idx])
            stations[name]['values'] = self.nc['zeta'][:, idx]
        self.__stations = stations

    def _init_netcdf_datetime(self):
        base_date = self.nc['time'].base_date.split('!')[0].strip(' ')
        for fmt in ('%Y-%m-%dT%H:%M', '%Y-%m-%d %H:%M'):
            try:
                base_date = datetime.strptime(base_date, fmt)
                break
            except ValueError:
                pass
        if isinstance(base_date, str):
            msg = f"Could not parse input date {base_date}. "
            msg += "Known formats are '%Y-%m-%dT%H:%M', '%Y-%m-%d %H:%M'."
            raise IOError(msg)
        self.__datetime = [base_date + timedelta(seconds=float(s))
                           for s in self.nc['time']]

    def _init_ascii(self):
        msg = "ASCII fort.61 files have not yet been implemented."
        raise NotImplementedError(msg)

    @property
    def path(self):
        return self._path

    @property
    def stations(self):
        return self._stations

    @property
    def datetime(self):
        return self._datetime

    @property
    def nc(self):
        return self._nc

    @property
    def _path(self):
        return self.__path

    @property
    def _stations(self):
        try:
            return self.__stations
        except AttributeError:
            if self.nc:
                self._init_netcdf_stations()
            else:
                self._init_ascii()
            return self.__stations

    @property
    def _nc(self):
        try:
            return self.__nc
        except AttributeError:
            try:
                nc = Dataset(self.path)
            except OSError:
                nc = False
            if nc:
                self._certify_netcdf_stations_file(nc)
            self.__nc = nc
            return self.__nc

    @property
    def _datetime(self):
        try:
            return self.__datetime
        except AttributeError:
            if self.nc:
                self._init_netcdf_datetime()
            else:
                self._init_ascii()
            return self.__datetime

    @_path.setter
    def _path(self, path):
        path = pathlib.Path(path)
        self.__path = path


# alias
Fort61 = ElevationStations
