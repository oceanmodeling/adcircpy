from abc import abstractmethod
from datetime import datetime, timedelta

from haversine import haversine
import numpy
from pyproj import CRS, Geod, Proj, Transformer


class WindForcing:
    def __init__(self, start_date: datetime, end_date: datetime, crs: CRS, nws: int = 20):
        self.start_date = start_date
        self.end_date = end_date
        self.NWS = nws

        if not isinstance(crs, CRS):
            crs = CRS.from_user_input(crs)
        self.__crs = crs

    @property
    def name(self):
        return self.df['name'].value_counts()[:].index.tolist()[0]

    @property
    def basin(self):
        return self.df['basin'].iloc[0]

    @property
    def storm_number(self):
        return self.df['storm_number'].iloc[0]

    @property
    def year(self):
        return self.df['datetime'].iloc[0].year

    @property
    def datetime(self):
        return self.df['datetime']

    @property
    def speed(self):
        return self.df['speed']

    @property
    def direction(self):
        return self.df['direction']

    @property
    def longitude(self):
        return self.df['longitude']

    @property
    def latitude(self):
        return self.df['latitude']

    @property
    def crs(self):
        return self.__crs

    @property
    def df(self):
        return self._df[(self._df['datetime'] >= self.start_date) &
                        (self._df['datetime'] <= self._file_end_date)]

    @property
    def start_date(self):
        return self._start_date

    @start_date.setter
    def start_date(self, start_date: datetime):
        self._start_date = start_date

    @property
    @abstractmethod
    def _start_date(self):
        raise NotImplementedError

    @_start_date.setter
    @abstractmethod
    def _start_date(self, start_date):
        raise NotImplementedError

    @property
    def end_date(self):
        return self._end_date

    @end_date.setter
    def end_date(self, end_date: datetime):
        self._end_date = end_date

    @property
    @abstractmethod
    def _end_date(self):
        raise NotImplementedError

    @_end_date.setter
    @abstractmethod
    def _end_date(self, end_date):
        raise NotImplementedError

    @property
    def _file_end_date(self):
        for date in numpy.unique(self._df['datetime']):
            if date >= self.end_date:
                return date

    def transform_to(self, crs: CRS):
        if self.crs != crs:
            # TODO implement reprojection
            pass

    @property
    def NWS(self):
        try:
            return self.__nws
        except AttributeError:
            return 20

    @NWS.setter
    def NWS(self, nws: int):
        self.__nws = int(nws)

    @property
    @abstractmethod
    def _df(self):
        raise NotImplementedError

    @property
    @abstractmethod
    def fort22(self):
        raise NotImplementedError

    @property
    @abstractmethod
    def WTIMINC(self):
        raise NotImplementedError

    @property
    @abstractmethod
    def BLADj(self):
        raise NotImplementedError

    @BLADj.setter
    @abstractmethod
    def BLADj(self, BLADj: float):
        raise NotImplementedError

    @property
    @abstractmethod
    def geofactor(self):
        raise NotImplementedError

    @geofactor.setter
    @abstractmethod
    def geofactor(self, geofactor: float):
        raise NotImplementedError

    @staticmethod
    def _compute_velocity(data: {}) -> {}:
        """ Output has units of meters per second. """

        merc = Proj('EPSG:3395')
        x, y = merc(data['longitude'], data['latitude'])
        t = data['datetime']
        unique_datetimes = numpy.unique(t)
        for i, _datetime in enumerate(unique_datetimes):
            indexes, = numpy.where(numpy.asarray(t) == _datetime)
            for idx in indexes:
                if indexes[-1] + 1 < len(t):
                    dx = haversine((y[idx], x[indexes[-1] + 1]), (y[idx], x[idx]), unit='nmi')
                    dy = haversine((y[indexes[-1] + 1], x[idx]), (y[idx], x[idx]), unit='nmi')
                    dt = ((t[indexes[-1] + 1] - t[idx]) / timedelta(hours=1))
                    vx = numpy.copysign(dx / dt, x[indexes[-1] + 1] - x[idx])
                    vy = numpy.copysign(dy / dt, y[indexes[-1] + 1] - y[idx])
                else:
                    dx = haversine((y[idx], x[indexes[0] - 1]), (y[idx], x[idx]), unit='nmi')
                    dy = haversine((y[indexes[0] - 1], x[idx]), (y[idx], x[idx]), unit='nmi')
                    dt = ((t[idx] - t[indexes[0] - 1]) / timedelta(hours=1))
                    vx = numpy.copysign(dx / dt, x[idx] - x[indexes[0] - 1])
                    vy = numpy.copysign(dy / dt, y[idx] - y[indexes[0] - 1])
                speed = numpy.sqrt(dx ** 2 + dy ** 2) / dt
                bearing = (360. + numpy.rad2deg(numpy.arctan2(vx, vy))) % 360
                data['speed'].append(int(numpy.around(speed, 0)))
                data['direction'].append(int(numpy.around(bearing, 0)))
        return data


def ellipsoidal_distance(point_a: (float, float), point_b: (float, float), crs_a: CRS,
                         crs_b: CRS = None) -> float:
    if crs_b is not None:
        transformer = Transformer.from_crs(crs_b, crs_a)
        point_b = transformer.transform(*point_b)
    datum_json = crs_a.datum.to_json_dict()
    ellipsoid = Geod(a=datum_json['ellipsoid']['semi_major_axis'],
                     rf=datum_json['ellipsoid']['inverse_flattening'])
    return ellipsoid.inv(point_a[0], point_a[1], point_b[0], point_b[1])[2]
