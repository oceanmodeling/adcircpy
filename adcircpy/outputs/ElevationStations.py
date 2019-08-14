from datetime import datetime, timedelta
from pathlib import Path
from netCDF4 import Dataset


class ElevationStations:

    def __init__(self, path):
        self.__path = str(Path(path))

    def __iter__(self):
        for key, item in self.stations.items():
            yield key, item

    def __repr__(self):
        print(self.stations)

    def __init_netcdf(self):
        def get_base_date():
            # print(self.nc['time'][:])
            base_date = self.nc['time'].base_date.split('!')[0].strip(' ')
            for fmt in ('%Y-%m-%dT%H:%M', '%Y-%m-%d %H:%M'):
                try:
                    return datetime.strptime(base_date, fmt)
                except ValueError:
                    pass
            raise ValueError('no valid date format found')
        msg = 'Input file {} '.format(self.path)
        msg + 'is not an ADCIRC stations output file (fort.61.nc).'
        assert 'station_name' in self.nc.variables \
            and 'zeta' in self.nc.variables, msg
        stations = dict()
        base_date = get_base_date()
        self.__datetime = [base_date + timedelta(seconds=float(s))
                           for s in self.nc['time']]
        for idx, name in enumerate(self.nc['station_name']):
            name = ''.join(
                [s.decode('UTF-8') for s in name]).strip(' ')
            stations[name] = dict()
            stations[name]['x'] = float(self.nc['x'][idx])
            stations[name]['y'] = float(self.nc['y'][idx])
            stations[name]['values'] = self.nc['zeta'][:, idx]
        self.__stations = stations

    def __init_ascii(self):
        # f = self.f
        raise NotImplementedError('Ascii files not yet implemented.')

    @property
    def stations(self):
        if not hasattr(self, "__stations"):
            try:
                self.__init_netcdf()
            except OSError:
                self.__init_ascii()
        return self.__stations

    @property
    def datetime(self):
        if not hasattr(self, "__datetime"):
            try:
                self.__init_netcdf()
            except OSError:
                self.__init_ascii()
        return self.__datetime

    @property
    def path(self):
        return self.__path

    @property
    def f(self):
        if not hasattr(self, "__f"):
            try:
                self.__f = open(self.path, 'r')
            except FileNotFoundError:
                raise
        return self.__f

    @property
    def nc(self):
        if not hasattr(self, "__nc"):
            try:
                self.__nc = Dataset(self.path)
            except OSError:
                raise
        return self.__nc


# alias
Fort61 = ElevationStations
