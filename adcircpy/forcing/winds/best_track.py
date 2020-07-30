from datetime import datetime, timedelta
import gzip
from io import BytesIO
import pathlib
import urllib.error
import urllib.request

# import utm
from haversine import haversine
from matplotlib import pyplot
from matplotlib.transforms import Bbox
import numpy
from pandas import DataFrame
from pyproj import CRS, Geod, Proj, Transformer
from shapely.geometry import Point, Polygon

from adcircpy.forcing.winds import atcf_id
from adcircpy.forcing.winds.base import WindForcing


# import os
# from pathlib import Path
# import zipfile
# from adcircpy.lib._get_cache_directory import _get_cache_directory


class BestTrackForcing(WindForcing):
    def __init__(self, storm_id: str, start_date: datetime = None, end_date: datetime = None,
                 crs: CRS = None):
        self._storm_id = storm_id
        super().__init__(start_date, end_date, crs)

    def clip_to_bbox(self, bbox: Bbox):
        """
        Important: bbox must be expressed in Mercator projection (EPSG:3395)
        """
        assert isinstance(bbox, Bbox), f"bbox must be a {Bbox} instance."
        bbox_pol = Polygon([
            [bbox.xmin, bbox.ymin],
            [bbox.xmax, bbox.ymin],
            [bbox.xmax, bbox.ymax],
            [bbox.xmin, bbox.ymax],
            [bbox.xmin, bbox.ymin]
        ])
        _switch = True
        unique_dates = numpy.unique(self._df['datetime'])
        for _datetime in unique_dates:
            records = self._df[self._df['datetime'] == _datetime]
            radii = records['radius_of_last_closed_isobar'].iloc[0]
            radii = 1852. * radii  # convert to meters
            merc = Proj("EPSG:3395")
            x, y = merc(
                records['longitude'].iloc[0],
                records['latitude'].iloc[0])
            p = Point(x, y)
            pol = p.buffer(radii)
            if _switch:
                if not pol.intersects(bbox_pol):
                    continue
                else:
                    self.start_date = records['datetime'].iloc[0]
                    _switch = False
                    continue
                    # self.start_date =
            else:
                if pol.intersects(bbox_pol):
                    continue
                else:
                    self.end_date = records['datetime'].iloc[0]
                    break

    def plot_trajectory(self, ax: pyplot.Axes = None, show: bool = False, color='k', **kwargs):
        kwargs.update({'color': color})
        if ax is None:
            fig = pyplot.figure()
            ax = fig.add_subplot(111)
        for i in range(len(self.speed)):
            # when dealing with nautical degrees, U is sine and V is cosine.
            U = self.speed.iloc[i] * numpy.sin(numpy.deg2rad(self.direction.iloc[i]))
            V = self.speed.iloc[i] * numpy.cos(numpy.deg2rad(self.direction.iloc[i]))
            ax.quiver(
                self.longitude.iloc[i], self.latitude.iloc[i], U, V, **kwargs)
            ax.annotate(
                self.df['datetime'].iloc[i],
                (self.longitude.iloc[i], self.latitude.iloc[i])
            )
        if show:
            ax.axis('scaled')
            pyplot.show()

    def write(self, path: str, overwrite: bool = False):
        path = pathlib.Path(path)
        if path.is_file() and not overwrite:
            raise Exception('Files exist, set overwrite=True to allow overwrite.')
        with open(path, 'w') as f:
            f.write(self.fort22)

    @property
    def storm_id(self):
        return self._storm_id

    @property
    def _storm_id(self):
        return f"{self.basin}{self.storm_number}{self.year}"

    @_storm_id.setter
    def _storm_id(self, storm_id):
        chars = 0
        for char in storm_id:
            if char.isdigit():
                chars += 1

        if chars == 4:
            _atcf_id = atcf_id(storm_id)
            if _atcf_id is None:
                raise Exception(f'No storm with id: {storm_id}')
            storm_id = _atcf_id

        url = f'ftp://ftp.nhc.noaa.gov/atcf/archive/{storm_id[4:]}/b{storm_id[0:2].lower()}{storm_id[2:]}.dat.gz'
        try:
            response = urllib.request.urlopen(url)
        except urllib.error.URLError:
            raise NameError(f'Did not find storm with id {storm_id} at url "{url}"')
        self.__atcf = BytesIO(response.read())

    @property
    def _start_date(self):
        return self.__start_date

    @_start_date.setter
    def _start_date(self, start_date):
        if start_date is not None:
            assert isinstance(start_date, datetime)
        else:
            start_date = self._df['datetime'].iloc[0]
        assert self._df['datetime'].iloc[0] <= start_date < self._df['datetime'].iloc[-1], \
            f"start_date must be {self._df['datetime'].iloc[0]} <= start_date ({start_date}) < " \
            f"{self._df['datetime'].iloc[-1]}"
        self.__start_date = start_date

    @property
    def _end_date(self):
        return self.__end_date

    @_end_date.setter
    def _end_date(self, end_date):
        if end_date is not None:
            assert isinstance(end_date, datetime)
        else:
            end_date = self._df['datetime'].iloc[-1]
        assert self._df['datetime'].iloc[0] < end_date <= self._df['datetime'].iloc[-1], \
            f"end_date must be {self._df['datetime'].iloc[0]} <= end_date ({end_date}) <= " \
            f"{self._df['datetime'].iloc[-1]}"
        assert end_date > self.start_date, \
            f"end_date ({end_date}) must be after start_date ({self.start_date})"
        self.__end_date = end_date

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
    def df(self):
        start_date_mask = self._df["datetime"] >= self.start_date
        end_date_mask = self._df["datetime"] <= self._file_end_date
        return self._df[start_date_mask & end_date_mask]

    @property
    def _df(self):
        # https://www.nrlmry.navy.mil/atcf_web/docs/database/new/abdeck.txt
        try:
            return self.__df
        except AttributeError:
            data = {
                "basin"                       : [],
                "storm_number"                : [],
                "datetime"                    : [],
                "record_type"                 : [],
                "latitude"                    : [],
                "longitude"                   : [],
                "max_sustained_wind_speed"    : [],
                "central_pressure"            : [],
                "development_level"           : [],
                "isotach"                     : [],
                "quadrant"                    : [],
                "radius_for_NEQ"              : [],
                "radius_for_SEQ"              : [],
                "radius_for_SWQ"              : [],
                "radius_for_NWQ"              : [],
                "background_pressure"         : [],
                "radius_of_last_closed_isobar": [],
                "radius_of_maximum_winds"     : [],
                "name"                        : [],
                "direction"                   : [],
                "speed"                       : []
            }
            for i, line in enumerate(gzip.GzipFile(fileobj=self.__atcf)):
                line = line.decode('UTF-8').split(',')
                data['basin'].append(line[0])
                data['storm_number'].append(line[1].strip(' '))
                _datetime = line[2].strip(' ')
                _minutes = line[3].strip(' ')
                if _minutes == '':
                    _minutes = '00'
                _datetime = _datetime + _minutes
                data['datetime'].append(datetime.strptime(_datetime, '%Y%m%d%H%M'))
                data['record_type'].append(line[4].strip(' '))
                if 'N' in line[6]:
                    _lat = float(line[6].strip('N ')) * 0.1
                elif 'S' in line:
                    _lat = float(line[6].strip('S ')) * -0.1
                data['latitude'].append(_lat)
                if 'E' in line[7]:
                    _lon = float(line[7].strip('E ')) * 0.1
                elif 'W' in line[7]:
                    _lon = float(line[7].strip('W ')) * -0.1
                data['longitude'].append(_lon)
                data['max_sustained_wind_speed'].append(float(line[8].strip(' ')))
                data['central_pressure'].append(float(line[9].strip(' ')))
                data['development_level'].append(line[10].strip(' '))
                data['isotach'].append(int(line[11].strip(' ')))
                data['quadrant'].append(line[12].strip(' '))
                data['radius_for_NEQ'].append(int(line[13].strip(' ')))
                data['radius_for_SEQ'].append(int(line[14].strip(' ')))
                data['radius_for_SWQ'].append(int(line[15].strip(' ')))
                data['radius_for_NWQ'].append(int(line[16].strip(' ')))
                if len(line) > 18:
                    data['background_pressure'].append(int(line[17].strip(' ')))
                    data['radius_of_last_closed_isobar'].append(int(line[18].strip(' ')))
                    data['radius_of_maximum_winds'].append(int(line[19].strip(' ')))
                    if len(line) > 23:
                        data['name'].append(line[27].strip(' '))
                    else:
                        data['name'].append('')
                else:
                    data['background_pressure'].append(data['background_pressure'][-1])
                    data['radius_of_last_closed_isobar'].append(
                        data['radius_of_last_closed_isobar'][-1])
                    data['radius_of_maximum_winds'].append(data['radius_of_maximum_winds'][-1])
                    data['name'].append('')
            data = self._compute_velocity(data)
            # data = self._transform_coordinates(data)
            self.__df = DataFrame(data=data)
            return self.__df

    @property
    def fort22(self):
        record_number = self._generate_record_numbers()
        fort22 = ''
        for i, (_, row) in enumerate(self.df.iterrows()):
            fort22 += f'{row["basin"]:<2},'
            fort22 += f'{row["storm_number"]:>3},'
            fort22 += f'{format(row["datetime"], "%Y%m%d%H"):>11},'
            fort22 += f'{"":3},'
            fort22 += f'{row["record_type"]:>5},'
            fort22 += f'{int((row["datetime"] - self.start_date) / timedelta(hours=1)):>4},'
            if row["latitude"] >= 0:
                fort22 += f'{int(row["latitude"] / 0.1):>4}N,'
            else:
                fort22 += f'{int(row["latitude"] / -0.1):>4}S,'
            if row["longitude"] >= 0:
                fort22 += f'{int(row["longitude"] / 0.1):>5}E,'
            else:
                fort22 += f'{int(row["longitude"] / -0.1):>5}W,'
            fort22 += f'{int(row["max_sustained_wind_speed"]):>4},'
            fort22 += f'{int(row["central_pressure"]):>5},'
            fort22 += f'{row["development_level"]:>3},'
            fort22 += f'{int(row["isotach"]):>4},'
            fort22 += f'{row["quadrant"]:>4},'
            fort22 += f'{int(row["radius_for_NEQ"]):>5},'
            fort22 += f'{int(row["radius_for_SEQ"]):>5},'
            fort22 += f'{int(row["radius_for_SWQ"]):>5},'
            fort22 += f'{int(row["radius_for_NWQ"]):>5},'
            if row["background_pressure"] is None:
                row["background_pressure"] = self.df["background_pressure"].iloc[i - 1]
            if row["background_pressure"] <= row["central_pressure"] < 1013:
                fort22 += f'{1013:>5},'
            elif row["background_pressure"] <= row["central_pressure"] >= 1013:
                fort22 += f'{int(row["central_pressure"] + 1):>5},'
            else:
                fort22 += f'{int(row["background_pressure"]):>5},'
            fort22 += f'{int(row["radius_of_last_closed_isobar"]):>5},'
            fort22 += f'{int(row["radius_of_maximum_winds"]):>4},'
            fort22 += f'{"":>5},'  # gust
            fort22 += f'{"":>4},'  # eye
            fort22 += f'{"":>4},'  # subregion
            fort22 += f'{"":>4},'  # maxseas
            fort22 += f'{"":>4},'  # initials
            fort22 += f'{row["direction"]:>3},'
            fort22 += f'{row["speed"]:>4},'
            fort22 += f'{row["name"]:^12},'
            # from this point forwards it's all aswip
            fort22 += f'{record_number[i]:>4},'
            fort22 += '\n'
        return fort22

    @property
    def WTIMINC(self):
        WTIMINC = self.start_date.strftime('%Y %m %d %H ')
        WTIMINC += f'{self.df["storm_number"].iloc[0]} '
        WTIMINC += f'{self.BLADj} '
        WTIMINC += f'{self.geofactor}'
        return WTIMINC

    @property
    def BLADj(self):
        try:
            return self.__BLADj
        except AttributeError:
            return 0.9

    @BLADj.setter
    def BLADj(self, BLADj: float):
        BLADj = float(BLADj)
        assert 0 <= BLADj <= 1
        self.__BLADj = BLADj

    @property
    def geofactor(self):
        try:
            return self.__geofactor
        except AttributeError:
            return 1

    @geofactor.setter
    def geofactor(self, geofactor: float):
        geofactor = float(geofactor)
        assert 0 <= geofactor <= 1
        self.__geofactor = geofactor

    def _generate_record_numbers(self):
        record_number = [1]
        for i in range(1, len(self.datetime)):
            if self.datetime.iloc[i] == self.datetime.iloc[i - 1]:
                record_number.append(record_number[-1])
            else:
                record_number.append(record_number[-1] + 1)
        return record_number

    @property
    def _file_end_date(self):
        unique_dates = numpy.unique(self._df['datetime'])
        for date in unique_dates:
            if date >= self.end_date:
                return date

    @staticmethod
    def _compute_velocity(data: {}):
        """
        Output has units of meters per second.
        """

        merc = Proj("EPSG:3395")
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
    if isinstance(point_a, Point):
        point_a = [*point_a.coords]
    if isinstance(point_b, Point):
        point_b = [*point_b.coords]
    if crs_b is not None:
        transformer = Transformer.from_crs(crs_b, crs_a)
        point_b = transformer.transform(*point_b)
    datum_json = crs_a.datum.to_json_dict()
    ellipsoid = Geod(a=datum_json['ellipsoid']['semi_major_axis'],
                     rf=datum_json['ellipsoid']['inverse_flattening'])
    return ellipsoid.inv(point_a[0], point_a[1], point_b[0], point_b[1])[2]
