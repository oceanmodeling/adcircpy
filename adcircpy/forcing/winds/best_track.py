from collections.abc import Collection
from datetime import datetime, timedelta
from functools import partial, wraps
import gzip
import io
from io import StringIO
import logging
import os
from os import PathLike
import pathlib
import time
from typing import Any, Union
import urllib.request
import zipfile

import appdirs
import geopandas
from haversine import haversine
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import numpy
import numpy as np
from pandas import DataFrame, read_csv
import pyproj
from pyproj import CRS, Proj
from shapely import ops
from shapely.geometry import Point, Polygon


from adcircpy.forcing.winds.base import WindForcing

logger = logging.getLogger(__name__)


def _dist(points1, points2) -> np.ndarray:
    R = 6373.0  # radius of earth
    lons1, lats1 = np.hsplit(points1, 2)
    lons2, lats2 = np.hsplit(points2, 2)
    lats1 = np.radians(lats1)
    lons1 = np.radians(lons1)
    lats2 = np.radians(lats2)
    lons2 = np.radians(lons2)

    dlons = lons2 - lons1
    dlats = lats2 - lats1

    a = np.sin(dlats / 2) ** 2 + np.cos(lats1) * np.cos(lats2) * np.sin(dlons / 2) ** 2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    return R * c


def _fetch_and_plot_coastline(_ax, show, filename='BestTrack.png'):
    save_dir = pathlib.Path(appdirs.user_data_dir('ne_coastline'))
    save_dir.mkdir(exist_ok=True, parents=True)

    file_path = save_dir / 'ne_110m_coastline.shp'
    zip_file_path = save_dir / 'ne_110m_coastline.zip'

    if not (save_dir / 'ne_110m_coastline.shp').exists():
        # download and save if not present
        url = 'http://naciscdn.org/naturalearth/'
        url += '110m/physical/ne_110m_coastline.zip'
        urllib.request.urlretrieve(url, zip_file_path)
        with zipfile.ZipFile(zip_file_path) as _zip:
            for name in _zip.namelist():
                data = _zip.read(name)
                outfile = os.path.join(save_dir, name)
                f = open(outfile, 'wb')
                f.write(data)
                f.close()
        os.remove(zip_file_path)

    open_shapefile = geopandas.read_file(file_path)

    open_shapefile.plot(ax=_ax)
    if show:
        plt.show()
    else:
        plt.savefig(filename)


class BestTrackForcing(WindForcing):
    def __init__(
        self,
        storm: Union[str, PathLike, DataFrame, io.BytesIO],
        nws: int = None,
        interval_seconds: int = None,
        start_date: datetime = None,
        end_date: datetime = None,
        dst_crs: CRS = None,
    ):
        if nws is None:
            nws = 20

        if interval_seconds is None:
            interval_seconds = 3600

        if isinstance(storm, DataFrame):
            self.__df = storm
        elif isinstance(storm, io.BytesIO):
            self.__atcf = storm
        elif isinstance(storm, (str, PathLike, pathlib.Path)):
            if os.path.exists(storm):
                self.__atcf = io.open(storm, 'rb')
            else:
                self._storm_id = storm

        assert nws in [8, 19, 20]
        self._start_date = start_date
        self._end_date = end_date
        self._dst_crs = dst_crs

        super().__init__(nws=nws, interval_seconds=interval_seconds)

    def summary(self, output: Union[str, os.PathLike] = None, overwrite: bool = False):
        min_storm_speed = np.min(self.speed)
        max_storm_speed = np.max(self.speed)
        track_length = self.track_length
        duration = self.duration
        min_central_pressure = np.min(self.central_pressure)
        max_wind_speed = np.max(self.df['max_sustained_wind_speed'])
        start_loc = (self.df['longitude'][0], self.df['latitude'][0])
        end_loc = (self.df['longitude'].iloc[-1], self.df['latitude'].iloc[-1])
        f = [
            f'Summary of storm: {self.storm_id}',
            f'min./max. track speed: {min_storm_speed} m/s, {max_storm_speed} m/s',
            f'min. central pressure: {min_central_pressure} hPa',
            f'max. wind speed: {max_wind_speed} kts',
            f'Starting at: {start_loc} and ended at: {end_loc}',
            f'Total track length: {track_length:.2f} km',
            f'Total track duration: {duration:.2f} days',
        ]
        summary = '\n'.join(f)
        if output is not None:
            if not isinstance(output, pathlib.Path):
                path = pathlib.Path(output)
            if path.exists() and overwrite is False:
                raise Exception('File exist, set overwrite=True to allow overwrite.')
            with open(path, 'w+') as fh:
                fh.write(summary)
        return summary

    def __str__(self):
        record_number = self._generate_record_numbers()
        fort22 = []
        for i, (_, row) in enumerate(self.df.iterrows()):
            line = []

            line.extend(
                [
                    f'{row["basin"]:<2}',
                    f'{row["storm_number"]:>3}',
                    f'{row["datetime"]:%Y%m%d%H}'.rjust(11),
                    f'{"":3}',
                    f'{row["record_type"]:>5}',
                    f'{convert_value((row["datetime"] - self.start_date) / timedelta(hours=1), to_type=int):>4}',
                ]
            )

            latitude = convert_value(row['latitude'] / 0.1, to_type=int, round_digits=1,)
            if latitude >= 0:
                line.append(f'{latitude:>4}N')
            else:
                line.append(f'{latitude * -.1:>4}S')

            longitude = convert_value(row['longitude'] / 0.1, to_type=int, round_digits=1,)
            if longitude >= 0:
                line.append(f'{longitude:>5}E')
            else:
                line.append(f'{longitude * -1:>5}W')

            line.extend(
                [
                    f'{convert_value(row["max_sustained_wind_speed"], to_type=int, round_digits=0):>4}',
                    f'{convert_value(row["central_pressure"], to_type=int, round_digits=0):>5}',
                    f'{row["development_level"]:>3}',
                    f'{convert_value(row["isotach"], to_type=int, round_digits=0):>4}',
                    f'{row["quadrant"]:>4}',
                    f'{convert_value(row["radius_for_NEQ"], to_type=int, round_digits=0):>5}',
                    f'{convert_value(row["radius_for_SEQ"], to_type=int, round_digits=0):>5}',
                    f'{convert_value(row["radius_for_SWQ"], to_type=int, round_digits=0):>5}',
                    f'{convert_value(row["radius_for_NWQ"], to_type=int, round_digits=0):>5}',
                ]
            )

            if row['background_pressure'] is None:
                row['background_pressure'] = self.df['background_pressure'].iloc[i - 1]
            if (
                row['background_pressure'] <= row['central_pressure']
                and 1013 > row['central_pressure']
            ):
                background_pressure = 1013
            elif (
                row['background_pressure'] <= row['central_pressure']
                and 1013 <= row['central_pressure']
            ):
                background_pressure = convert_value(
                    row['central_pressure'] + 1, to_type=int, round_digits=0,
                )
            else:
                background_pressure = convert_value(
                    row['background_pressure'], to_type=int, round_digits=0,
                )
            line.append(f'{background_pressure:>5}')

            line.extend(
                [
                    f'{convert_value(row["radius_of_last_closed_isobar"], to_type=int, round_digits=0):>5}',
                    f'{convert_value(row["radius_of_maximum_winds"], to_type=int, round_digits=0):>4}',
                    f'{"":>5}',  # gust
                    f'{"":>4}',  # eye
                    f'{"":>4}',  # subregion
                    f'{"":>4}',  # maxseas
                    f'{"":>4}',  # initials
                    f'{row["direction"]:>3}',
                    f'{row["speed"]:>4}',
                    f'{row["name"]:^12}',
                ]
            )

            # from this point forwards it's all aswip
            line.append(f'{record_number[i]:>4}')

            fort22.append(','.join(line))

        return '\n'.join(fort22)

    def write(self, path: PathLike, overwrite: bool = False):
        if not isinstance(path, pathlib.Path):
            path = pathlib.Path(path)
        if path.exists() and overwrite is False:
            raise Exception('File exist, set overwrite=True to allow overwrite.')
        with open(path, 'w') as f:
            f.write(str(self))

    @property
    def NWS(self) -> int:
        try:
            return self.__NWS
        except AttributeError:
            return 20

    @NWS.setter
    def NWS(self, NWS: int):
        assert NWS in [8, 19, 20]
        self.__NWS = int(NWS)

    @property
    def WTIMINC(self) -> str:
        return (
            f'{self.start_date:%Y %m %d %H} '
            f'{self.df["storm_number"].iloc[0]} '
            f'{self.BLADj} '
            f'{self.geofactor}'
        )

    @property
    def BLADj(self) -> float:
        try:
            return self.__BLADj
        except AttributeError:
            return 0.9

    @BLADj.setter
    def BLADj(self, BLADj: float):
        BLADj = float(BLADj)
        assert BLADj >= 0 and BLADj <= 1
        self.__BLADj = BLADj

    @property
    def geofactor(self) -> float:
        try:
            return self.__geofactor
        except AttributeError:
            return 1

    @geofactor.setter
    def geofactor(self, geofactor: float):
        geofactor = float(geofactor)
        assert geofactor >= 0 and geofactor <= 1
        self.__geofactor = geofactor

    @property
    def storm_id(self) -> str:
        return self._storm_id

    @property
    def _storm_id(self) -> str:
        return f'{self.basin}{self.storm_number}{self.year}'

    @_storm_id.setter
    def _storm_id(self, storm_id: str):
        if storm_id is not None:
            chars = 0
            for char in storm_id:
                if char.isdigit():
                    chars += 1

            if chars == 4:

                _atcf_id = atcf_id(storm_id)
                if _atcf_id is None:
                    msg = f'No storm with id: {storm_id}'
                    raise Exception(msg)
                storm_id = _atcf_id

            url = f'ftp://ftp.nhc.noaa.gov/atcf/archive/{storm_id[4:]}/b{storm_id[0:2].lower()}{storm_id[2:]}.dat.gz'

            try:
                logger.info(f'Downloading storm data from {url}')
                response = urllib.request.urlopen(url)
            except urllib.error.URLError as e:
                if '550' in e.reason:
                    raise NameError(
                        f'Did not find storm with id {storm_id}. '
                        + f'Submitted URL was {url}.'
                    )
                else:

                    @retry(urllib.error.URLError, tries=4, delay=3, backoff=2)
                    def make_request():
                        logger.info(f'Downloading storm data from {url} failed, retrying...')
                        return urllib.request.urlopen(url)

                    response = make_request()

            self.__atcf = io.BytesIO(response.read())

    @property
    def track_length(self) -> float:
        if not hasattr(self, '_track_length'):
            lons, lats = self.df['longitude'], self.df['latitude']
            prev = np.asarray([lons.iloc[:-1], lats.iloc[:-1]]).T
            curr = np.asarray([lons.iloc[1:], lats.iloc[1:]]).T
            distances = _dist(prev, curr)
            self._track_length = np.sum(distances)
        return self._track_length

    @property
    def start_date(self) -> datetime:
        return self._start_date

    @start_date.setter
    def start_date(self, start_date: datetime):
        self._start_date = start_date

    @property
    def _start_date(self) -> datetime:
        return self.__start_date

    @_start_date.setter
    def _start_date(self, start_date: datetime):
        if start_date is not None:
            assert isinstance(start_date, datetime)
        else:
            start_date = self._df['datetime'].iloc[0]
        msg = f"start_date must be >= {self._df['datetime'].iloc[0]} "
        msg += f"and <{self._df['datetime'].iloc[-1]}"
        assert (
            start_date >= self._df['datetime'].iloc[0]
            and start_date < self._df['datetime'].iloc[-1]
        ), msg
        self.__start_date = start_date

    @property
    def end_date(self) -> datetime:
        return self._end_date

    @end_date.setter
    def end_date(self, end_date: datetime):
        self._end_date = end_date

    @property
    def _end_date(self) -> datetime:
        return self.__end_date

    @_end_date.setter
    def _end_date(self, end_date: datetime):
        if end_date is not None:
            assert isinstance(end_date, datetime)
        else:
            end_date = self._df['datetime'].iloc[-1]
        msg = f"end_date must be >= {self._df['datetime'].iloc[0]} "
        msg += f"and <= {self._df['datetime'].iloc[-1]}. "
        msg += f'The given end_date was {end_date}.'
        assert (
            end_date > self._df['datetime'].iloc[0]
            and end_date <= self._df['datetime'].iloc[-1]
        ), msg
        msg = 'end_date must be larger than start_date.\n'
        msg += f'start_date is {self.start_date} and end_date is {end_date}.'
        assert end_date > self.start_date, msg
        self.__end_date = end_date

    @property
    def duration(self) -> float:
        d = (self.__end_date - self.__start_date).days
        return d

    @property
    def name(self) -> str:
        return self.df['name'].value_counts()[:].index.tolist()[0]

    @property
    def basin(self) -> str:
        return self.df['basin'].iloc[0]

    @property
    def storm_number(self) -> str:
        return self.df['storm_number'].iloc[0]

    @property
    def year(self) -> int:
        return self.df['datetime'].iloc[0].year

    @property
    def datetime(self):
        return self.df['datetime']

    @property
    def central_pressure(self):
        return self.df['central_pressure']

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
        start_date_mask = self._df['datetime'] >= self.start_date
        end_date_mask = self._df['datetime'] <= self._file_end_date
        return self._df[start_date_mask & end_date_mask]

    @property
    def _df(self):
        # https://www.nrlmry.navy.mil/atcf_web/docs/database/new/abdeck.txt
        try:
            return self.__df
        except AttributeError:
            data = {
                'basin': list(),
                'storm_number': list(),
                'datetime': list(),
                'record_type': list(),
                'latitude': list(),
                'longitude': list(),
                'max_sustained_wind_speed': list(),
                'central_pressure': list(),
                'development_level': list(),
                'isotach': list(),
                'quadrant': list(),
                'radius_for_NEQ': list(),
                'radius_for_SEQ': list(),
                'radius_for_SWQ': list(),
                'radius_for_NWQ': list(),
                'background_pressure': list(),
                'radius_of_last_closed_isobar': list(),
                'radius_of_maximum_winds': list(),
                'name': list(),
                'direction': list(),
                'speed': list(),
            }
            if isinstance(self.__atcf, io.BytesIO):
                fd = gzip.GzipFile(fileobj=self.__atcf)
            else:
                fd = self.__atcf

            for i, line in enumerate(fd):
                line = line.decode('UTF-8').split(',')
                data['basin'].append(line[0])
                data['storm_number'].append(line[1].strip(' '))
                _datetime = line[2].strip(' ')
                _minutes = line[3].strip(' ')
                if _minutes == "":
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
                try:
                    data['isotach'].append(int(line[11].strip(' ')))
                except ValueError:
                    raise Exception(
                        'Error: No radial wind information for this storm; '
                        'parametric wind model cannot be built.'
                    )
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
                        data['name'].append("")
                else:
                    data['background_pressure'].append(data['background_pressure'][-1])
                    data['radius_of_last_closed_isobar'].append(
                        data['radius_of_last_closed_isobar'][-1]
                    )
                    data['radius_of_maximum_winds'].append(data['radius_of_maximum_winds'][-1])
                    data['name'].append("")
            data = self._compute_velocity(data)
            # data = self._transform_coordinates(data)
            self.__df = DataFrame(data=data)
            return self.__df

    @_df.setter
    def _df(self, data):
        msg = f'data must be a pandas {DataFrame} instance.'
        assert isinstance(data, DataFrame), msg
        self.__df = DataFrame(data=data)

    def clip_to_bbox(self, bbox, bbox_crs):
        msg = f'bbox must be a {Bbox} instance.'
        assert isinstance(bbox, Bbox), msg
        bbox_pol = Polygon(
            [
                [bbox.xmin, bbox.ymin],
                [bbox.xmax, bbox.ymin],
                [bbox.xmax, bbox.ymax],
                [bbox.xmin, bbox.ymax],
                [bbox.xmin, bbox.ymin],
            ]
        )
        _switch = True
        unique_dates = np.unique(self._df['datetime'])
        _found_start_date = False
        for _datetime in unique_dates:
            records = self._df[self._df['datetime'] == _datetime]
            radii = records['radius_of_last_closed_isobar'].iloc[0]
            radii = 1852.0 * radii  # convert to meters
            lon = records['longitude'].iloc[0]
            lat = records['latitude'].iloc[0]
            pol = get_circle_of_radius(lon, lat, radii)
            if _switch is True:
                if not pol.intersects(bbox_pol):
                    continue
                else:
                    self.start_date = records['datetime'].iloc[0]
                    _found_start_date = True
                    _switch = False
                    continue

            else:
                if pol.intersects(bbox_pol):
                    continue
                else:
                    self.end_date = records['datetime'].iloc[0]
                    break

        if _found_start_date is False:
            raise Exception(f'No data within mesh bounding box for storm {self.storm_id}.')

    def plot_track(self, axes=None, show=False, color='k', annotate=True, **kwargs):
        kwargs.update({'color': color})
        if axes is None:
            fig = plt.figure()
            axes = fig.add_subplot(111)
        for i in range(len(self.speed)):
            # when dealing with nautical degrees, U is sine and V is cosine.
            U = self.speed.iloc[i] * np.sin(np.deg2rad(self.direction.iloc[i]))
            V = self.speed.iloc[i] * np.cos(np.deg2rad(self.direction.iloc[i]))
            axes.quiver(self.longitude.iloc[i], self.latitude.iloc[i], U, V, **kwargs)
            if annotate:
                if i % 6 == 0:
                    axes.annotate(
                        self.df['datetime'].iloc[i],
                        (self.longitude.iloc[i], self.latitude.iloc[i]),
                    )
        if show:
            axes.axis('scaled')
        _fetch_and_plot_coastline(axes, show)

    def _generate_record_numbers(self):
        record_number = [1]
        for i in range(1, len(self.datetime)):
            if self.datetime.iloc[i] == self.datetime.iloc[i - 1]:
                record_number.append(record_number[-1])
            else:
                record_number.append(record_number[-1] + 1)
        return record_number

    def transform_to(self, crs):
        pass

    @property
    def _file_end_date(self):
        unique_dates = np.unique(self._df['datetime'])
        for date in unique_dates:
            if date >= np.datetime64(self.end_date):
                return date

    @staticmethod
    def _compute_velocity(data):
        """
        Output has units of meters per second.
        """
        merc = Proj('EPSG:3395')
        x, y = merc(data['longitude'], data['latitude'])
        unique_datetimes = np.unique(data['datetime'])
        for i, _datetime in enumerate(unique_datetimes):
            (indexes,) = np.where(np.asarray(data['datetime']) == _datetime)
            for idx in indexes:
                if indexes[-1] + 1 < len(data['datetime']):
                    dt = (
                        data['datetime'][indexes[-1] + 1] - data['datetime'][idx]
                    ).total_seconds() / (60.0 * 60.0)
                    dx = haversine(
                        (data['latitude'][idx], data['longitude'][indexes[-1] + 1]),
                        (data['latitude'][idx], data['longitude'][idx]),
                        unit='nmi',
                    )
                    dy = haversine(
                        (data['latitude'][indexes[-1] + 1], data['longitude'][idx]),
                        (data['latitude'][idx], data['longitude'][idx]),
                        unit='nmi',
                    )
                    vx = np.copysign(
                        dx / dt, data['longitude'][indexes[-1] + 1] - data['longitude'][idx],
                    )
                    vy = np.copysign(
                        dy / dt, data['latitude'][indexes[-1] + 1] - data['latitude'][idx],
                    )
                else:
                    dt = (
                        data['datetime'][idx] - data['datetime'][indexes[0] - 1]
                    ).total_seconds() / (60.0 * 60.0)
                    dx = haversine(
                        (data['latitude'][idx], data['longitude'][indexes[0] - 1]),
                        (data['latitude'][idx], data['longitude'][idx]),
                        unit='nmi',
                    )
                    dy = haversine(
                        (data['latitude'][indexes[0] - 1], data['longitude'][idx]),
                        (data['latitude'][idx], data['longitude'][idx]),
                        unit='nmi',
                    )
                    vx = np.copysign(
                        dx / dt, data['longitude'][idx] - data['longitude'][indexes[0] - 1],
                    )
                    vy = np.copysign(
                        dy / dt, data['latitude'][idx] - data['latitude'][indexes[0] - 1],
                    )
                speed = np.sqrt(dx ** 2 + dy ** 2) / dt
                bearing = (360.0 + np.rad2deg(np.arctan2(vx, vy))) % 360
                data['speed'].append(int(np.around(speed, 0)))
                data['direction'].append(int(np.around(bearing, 0)))
        return data

    @classmethod
    def from_fort22(
        cls,
        fort22: PathLike,
        nws: int = None,
        start_date: datetime = None,
        end_date: datetime = None,
    ) -> 'BestTrackForcing':
        if nws is None:
            nws = 20

        data = read_atcf(fort22)

        storm_id = f'{data["name"][0]}{data["datetime"][0]:%Y}'

        if start_date is None:
            start_date = min(data['datetime'])
        if end_date is None:
            end_date = max(data['datetime'])

        instance = cls(storm=storm_id, nws=nws, start_date=start_date, end_date=end_date,)

        instance.__df = data

        return instance

    @classmethod
    def from_atcf_file(
        cls, atcf: PathLike, nws: int, start_date: datetime = None, end_date: datetime = None,
    ) -> 'BestTrackForcing':
        return cls(storm=atcf, nws=nws, start_date=start_date, end_date=end_date,)


def convert_value(value: Any, to_type: type, round_digits: int = None) -> Any:
    if value is not None and value != "":
        if not isinstance(value, str) and numpy.isinf(value):
            return value
        if round_digits is not None and issubclass(to_type, (int, float)):
            if isinstance(value, str):
                value = float(value)
            value = round(value, round_digits)
        value = to_type(value)
    return value


def retry(ExceptionToCheck, tries=4, delay=3, backoff=2, logger=None):
    """Retry calling the decorated function using an exponential backoff.

    http://www.saltycrane.com/blog/2009/11/trying-out-retry-decorator-python/
    original from: http://wiki.python.org/moin/PythonDecoratorLibrary#Retry

    :param ExceptionToCheck: the exception to check. may be a tuple of
        exceptions to check
    :type ExceptionToCheck: Exception or tuple
    :param tries: number of times to try (not retry) before giving up
    :type tries: int
    :param delay: initial delay between retries in seconds
    :type delay: int
    :param backoff: backoff multiplier e.g. value of 2 will double the delay
        each retry
    :type backoff: int
    :param logger: logger to use. If None, print
    :type logger: logging.Logger instance
    """

    def deco_retry(f):
        @wraps(f)
        def f_retry(*args, **kwargs):
            mtries, mdelay = tries, delay
            while mtries > 1:
                try:
                    return f(*args, **kwargs)
                except ExceptionToCheck as e:
                    msg = '%s, Retrying in %d seconds...' % (str(e), mdelay)
                    if logger:
                        logger.warning(msg)
                    # else:
                    #     print(msg)
                    time.sleep(mdelay)
                    mtries -= 1
                    mdelay *= backoff
            return f(*args, **kwargs)

        return f_retry  # true decorator

    return deco_retry


def atcf_id(storm_id):
    url = 'ftp://ftp.nhc.noaa.gov/atcf/archive/storm.table'

    @retry(urllib.error.URLError, tries=4, delay=3, backoff=2)
    def request_url():
        logger.info(f'Querying storm name from: {url}')
        return urllib.request.urlopen(url)

    res = request_url()
    df = read_csv(StringIO("".join([_.decode('utf-8') for _ in res])), header=None,)
    name = storm_id[:-4]
    year = storm_id[-4:]
    entry = df.loc[(df[0] == name.upper().rjust(10)) & (df[8] == int(year))]
    if len(entry) == 0:
        return None
    else:
        return entry[20].tolist()[0].strip()


def read_atcf(track: PathLike) -> DataFrame:
    with open(track) as track_file:
        track = track_file.readlines()

    data = {
        'basin': [],
        'storm_number': [],
        'datetime': [],
        'record_type': [],
        'latitude': [],
        'longitude': [],
        'max_sustained_wind_speed': [],
        'central_pressure': [],
        'development_level': [],
        'isotach': [],
        'quadrant': [],
        'radius_for_NEQ': [],
        'radius_for_SEQ': [],
        'radius_for_SWQ': [],
        'radius_for_NWQ': [],
        'background_pressure': [],
        'radius_of_last_closed_isobar': [],
        'radius_of_maximum_winds': [],
        'name': [],
        'direction': [],
        'speed': [],
    }

    for index, row in enumerate(track):
        row = [value.strip() for value in row.split(',')]

        row_data = {key: None for key in data}

        row_data['basin'] = row[0]
        row_data['storm_number'] = row[1]
        row_data['datetime'] = datetime.strptime(row[2], '%Y%m%d%H')
        row_data['record_type'] = row[4]

        latitude = row[6]
        if 'N' in latitude:
            latitude = float(latitude[:-1]) * 0.1
        elif 'S' in latitude:
            latitude = float(latitude[:-1]) * -0.1
        row_data['latitude'] = latitude

        longitude = row[7]
        if 'E' in longitude:
            longitude = float(longitude[:-1]) * 0.1
        elif 'W' in longitude:
            longitude = float(longitude[:-1]) * -0.1
        row_data['longitude'] = longitude

        row_data['max_sustained_wind_speed'] = convert_value(
            row[8], to_type=int, round_digits=0,
        )
        row_data['central_pressure'] = convert_value(row[9], to_type=int, round_digits=0,)
        row_data['development_level'] = row[10]
        row_data['isotach'] = convert_value(row[11], to_type=int, round_digits=0,)
        row_data['quadrant'] = row[12]
        row_data['radius_for_NEQ'] = convert_value(row[13], to_type=int, round_digits=0,)
        row_data['radius_for_SEQ'] = convert_value(row[14], to_type=int, round_digits=0,)
        row_data['radius_for_SWQ'] = convert_value(row[15], to_type=int, round_digits=0,)
        row_data['radius_for_NWQ'] = convert_value(row[16], to_type=int, round_digits=0,)
        row_data['background_pressure'] = convert_value(row[17], to_type=int, round_digits=0,)
        row_data['radius_of_last_closed_isobar'] = convert_value(
            row[18], to_type=int, round_digits=0,
        )
        row_data['radius_of_maximum_winds'] = convert_value(
            row[19], to_type=int, round_digits=0,
        )
        row_data['direction'] = row[25]
        row_data['speed'] = row[26]
        row_data['name'] = row[27]

        for key, value in row_data.items():
            if isinstance(data[key], Collection):
                data[key].append(value)
            elif data[key] is None:
                data[key] = value

    return DataFrame(data=data)


def get_circle_of_radius(lon, lat, radius):

    local_azimuthal_projection = '+proj=aeqd +R=6371000 +units=m ' f'+lat_0={lat} +lon_0={lon}'
    wgs84_to_aeqd = partial(
        pyproj.transform,
        pyproj.Proj('+proj=longlat +datum=WGS84 +no_defs'),
        pyproj.Proj(local_azimuthal_projection),
    )
    aeqd_to_wgs84 = partial(
        pyproj.transform,
        pyproj.Proj(local_azimuthal_projection),
        pyproj.Proj('+proj=longlat +datum=WGS84 +no_defs'),
    )

    center = Point(float(lon), float(lat))
    point_transformed = ops.transform(wgs84_to_aeqd, center)

    return ops.transform(aeqd_to_wgs84, point_transformed.buffer(radius))
