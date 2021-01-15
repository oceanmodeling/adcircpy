from datetime import datetime, timedelta
from functools import wraps
import gzip
import io
from io import StringIO
from os import PathLike
import pathlib
import time
import urllib.request

from haversine import haversine
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import numpy as np
from pandas import DataFrame, read_csv
from pyproj import Proj
from shapely.geometry import Point, Polygon

from adcircpy.forcing.winds.base import WindForcing


class BestTrackForcing(WindForcing):
    def __init__(self, storm_id, nws: int = 20, interval_seconds: int = 3600,
                 start_date=None, end_date=None, dst_crs=None):
        assert nws in [19, 20]
        super().__init__(nws, interval_seconds)
        self._storm_id = storm_id
        self._start_date = start_date
        self._end_date = end_date
        self._dst_crs = dst_crs

    def write(self, path: PathLike, overwrite: bool = False):
        if not isinstance(path, pathlib.Path):
            path = pathlib.Path(path)
        if path.exists() or overwrite:
            with open(path, 'w') as f:
                f.write(self.fort22)
        else:
            raise Exception('Files exist, '
                            'set overwrite=True to allow overwrite.')

    @property
    def fort22(self):
        record_number = self._generate_record_numbers()
        fort22 = ''
        for i, (_, row) in enumerate(self.df.iterrows()):
            fort22 += f'{row["basin"]:<2},{row["storm_number"]:>3},' \
                      f'{row["datetime"]:%Y%m%d%H>11},' \
                      f'{"":3},{row["record_type"]:>5},' \
                      f'{int((row["datetime"] - self.start_date) / timedelta(hours=1)):>4},'
            if row["latitude"] >= 0:
                fort22 += f'{int(row["latitude"] / .1):>4}N,'
            else:
                fort22 += f'{int(row["latitude"] / -.1):>4}S,'
            if row["longitude"] >= 0:
                fort22 += f'{int(row["longitude"] / .1):>5}E,'
            else:
                fort22 += f'{int(row["longitude"] / -.1):>5}W,'
            fort22 += f'{int(row["max_sustained_wind_speed"]):>4},' \
                      f'{int(row["central_pressure"]):>5},' \
                      f'{row["development_level"]:>3},' \
                      f'{int(row["isotach"]):>4},' \
                      f'{row["quadrant"]:>4},' \
                      f'{int(row["radius_for_NEQ"]):>5},' \
                      f'{int(row["radius_for_SEQ"]):>5},' \
                      f'{int(row["radius_for_SWQ"]):>5},' \
                      f'{int(row["radius_for_NWQ"]):>5},'
            if row["background_pressure"] is None:
                row["background_pressure"] = \
                    self.df["background_pressure"].iloc[i - 1]
            if (row["background_pressure"] <= row["central_pressure"]
                    and 1013 > row["central_pressure"]):
                fort22 += f'{1013:>5},'
            elif (row["background_pressure"] <= row["central_pressure"]
                  and 1013 <= row["central_pressure"]):
                fort22 += f'{int(row["central_pressure"] + 1):>5},'
            else:
                fort22 += f'{int(row["background_pressure"]):>5},'
            fort22 += f'{int(row["radius_of_last_closed_isobar"]):>5},' \
                      f'{int(row["radius_of_maximum_winds"]):>4},'
            fort22 += f'{"":>5},'  # gust
            fort22 += f'{"":>4},'  # eye
            fort22 += f'{"":>4},'  # subregion
            fort22 += f'{"":>4},'  # maxseas
            fort22 += f'{"":>4},'  # initials
            fort22 += f'{row["direction"]:>3},' \
                      f'{row["speed"]:>4},' \
                      f'{row["name"]:^12},'
            # from this point forwards it's all aswip
            fort22 += f'{record_number[i]:>4},' \
                      f'\n'
        return fort22

    @property
    def NWS(self) -> int:
        try:
            return self.__NWS
        except AttributeError:
            return 20

    @NWS.setter
    def NWS(self, NWS: int):
        assert NWS in [19, 20]
        self.__NWS = int(NWS)

    @property
    def WTIMINC(self) -> str:
        return f'{self.start_date:%Y %m %d %H} ' \
               f'{self.df["storm_number"].iloc[0]} ' \
               f'{self.BLADj} ' \
               f'{self.geofactor}'

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
                msg = f'No storm with id: {storm_id}'
                raise Exception(msg)
            storm_id = _atcf_id

        url = 'ftp://ftp.nhc.noaa.gov/atcf/archive/'
        url += storm_id[4:]
        url += '/b'
        url += storm_id[0:2].lower()
        url += storm_id[2:]
        url += '.dat.gz'

        try:
            response = urllib.request.urlopen(url)
        except urllib.error.URLError as e:
            if '550' in e.reason:
                raise NameError(
                        f'Did not find storm with id {storm_id}. '
                        + f'Submitted URL was {url}.')
            else:
                @retry(urllib.error.URLError, tries=4, delay=3, backoff=2)
                def make_request():
                    return urllib.request.urlopen(url)

                response = make_request()

        self.__atcf = io.BytesIO(response.read())

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
        assert start_date >= self._df['datetime'].iloc[0] \
               and start_date < self._df['datetime'].iloc[-1], msg
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
        msg += f"The given end_date was {end_date}."
        assert end_date > self._df['datetime'].iloc[0] \
               and end_date <= self._df['datetime'].iloc[-1], msg
        msg = "end_date must be larger than start_date.\n"
        msg += f"start_date is {self.start_date} and end_date is {end_date}."
        assert end_date > self.start_date, msg
        self.__end_date = end_date

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
                "basin": list(),
                "storm_number": list(),
                "datetime": list(),
                "record_type": list(),
                "latitude": list(),
                "longitude": list(),
                "max_sustained_wind_speed": list(),
                "central_pressure": list(),
                "development_level": list(),
                "isotach": list(),
                "quadrant": list(),
                "radius_for_NEQ": list(),
                "radius_for_SEQ": list(),
                "radius_for_SWQ": list(),
                "radius_for_NWQ": list(),
                "background_pressure": list(),
                "radius_of_last_closed_isobar": list(),
                "radius_of_maximum_winds": list(),
                "name": list(),
                "direction": list(),
                "speed": list()
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
                data['datetime'].append(
                        datetime.strptime(_datetime, '%Y%m%d%H%M'))
                data['record_type'].append(line[4].strip(' '))
                if 'N' in line[6]:
                    _lat = float(line[6].strip('N ')) * .1
                elif 'S' in line:
                    _lat = float(line[6].strip('S ')) * -.1
                data['latitude'].append(_lat)
                if 'E' in line[7]:
                    _lon = float(line[7].strip('E ')) * .1
                elif 'W' in line[7]:
                    _lon = float(line[7].strip('W ')) * -.1
                data['longitude'].append(_lon)
                data['max_sustained_wind_speed'].append(
                        float(line[8].strip(' ')))
                data['central_pressure'].append(
                        float(line[9].strip(' ')))
                data['development_level'].append(
                        line[10].strip(' '))
                try:
                    data['isotach'].append(int(line[11].strip(' ')))
                except ValueError:
                    raise Exception(
                            'Error: No radial wind information for this storm; '
                            'parametric wind model cannot be built.')
                data['quadrant'].append(line[12].strip(' '))
                data['radius_for_NEQ'].append(
                        int(line[13].strip(' ')))
                data['radius_for_SEQ'].append(
                        int(line[14].strip(' ')))
                data['radius_for_SWQ'].append(
                        int(line[15].strip(' ')))
                data['radius_for_NWQ'].append(
                        int(line[16].strip(' ')))
                if len(line) > 18:
                    data['background_pressure'].append(
                            int(line[17].strip(' ')))
                    data['radius_of_last_closed_isobar'].append(
                            int(line[18].strip(' ')))
                    data['radius_of_maximum_winds'].append(
                            int(line[19].strip(' ')))
                    if len(line) > 23:
                        data['name'].append(line[27].strip(' '))
                    else:
                        data['name'].append('')
                else:
                    data['background_pressure'].append(
                            data['background_pressure'][-1])
                    data['radius_of_last_closed_isobar'].append(
                            data['radius_of_last_closed_isobar'][-1])
                    data['radius_of_maximum_winds'].append(
                            data['radius_of_maximum_winds'][-1])
                    data['name'].append('')
            data = self._compute_velocity(data)
            # data = self._transform_coordinates(data)
            self.__df = DataFrame(data=data)
            return self.__df

    def clip_to_bbox(self, bbox):
        """
        Important: bbox must be expressed in Mercator projection (EPSG:3395)
        """
        msg = f"bbox must be a {Bbox} instance."
        assert isinstance(bbox, Bbox), msg
        bbox_pol = Polygon(
                [[bbox.xmin, bbox.ymin],
                 [bbox.xmax, bbox.ymin],
                 [bbox.xmax, bbox.ymax],
                 [bbox.xmin, bbox.ymax],
                 [bbox.xmin, bbox.ymin]
                 ])
        _switch = True
        unique_dates = np.unique(self._df['datetime'])
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

    def plot_trajectory(self, ax=None, show=False, color='k', **kwargs):
        kwargs.update({'color': color})
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        for i in range(len(self.speed)):
            # when dealing with nautical degrees, U is sine and V is cosine.
            U = self.speed.iloc[i] * np.sin(np.deg2rad(self.direction.iloc[i]))
            V = self.speed.iloc[i] * np.cos(np.deg2rad(self.direction.iloc[i]))
            ax.quiver(
                    self.longitude.iloc[i], self.latitude.iloc[i], U, V,
                    **kwargs)
            ax.annotate(
                    self.df['datetime'].iloc[i],
                    (self.longitude.iloc[i], self.latitude.iloc[i])
            )
        if show:
            ax.axis('scaled')
            plt.show()

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
            if date >= self.end_date:
                return date

    @staticmethod
    def _compute_velocity(data):
        """
        Output has units of meters per second.
        """
        merc = Proj("EPSG:3395")
        x, y = merc(data['longitude'], data['latitude'])
        unique_datetimes = np.unique(data['datetime'])
        for i, _datetime in enumerate(unique_datetimes):
            indexes, = np.where(
                    np.asarray(data['datetime']) == _datetime)
            for idx in indexes:
                if indexes[-1] + 1 < len(data['datetime']):
                    dt = ((data['datetime'][indexes[-1] + 1]
                           - data['datetime'][idx])
                          .total_seconds() / (60. * 60.))
                    dx = haversine(
                            (data['latitude'][idx],
                             data['longitude'][indexes[-1] + 1]),
                            (data['latitude'][idx],
                             data['longitude'][idx]), unit='nmi')
                    dy = haversine(
                            (data['latitude'][indexes[-1] + 1],
                             data['longitude'][idx]),
                            (data['latitude'][idx],
                             data['longitude'][idx]), unit='nmi')
                    vx = np.copysign(
                            dx / dt,
                            data['longitude'][indexes[-1] + 1]
                            - data['longitude'][idx])
                    vy = np.copysign(
                            dy / dt,
                            data['latitude'][indexes[-1] + 1]
                            - data['latitude'][idx])
                else:
                    dt = ((data['datetime'][idx]
                           - data['datetime'][indexes[0] - 1])
                          .total_seconds() / (60. * 60.))
                    dx = haversine(
                            (data['latitude'][idx],
                             data['longitude'][indexes[0] - 1]),
                            (data['latitude'][idx],
                             data['longitude'][idx]), unit='nmi')
                    dy = haversine(
                            (data['latitude'][indexes[0] - 1],
                             data['longitude'][idx]),
                            (data['latitude'][idx],
                             data['longitude'][idx]), unit='nmi')
                    vx = np.copysign(
                            dx / dt,
                            data['longitude'][idx]
                            - data['longitude'][indexes[0] - 1])
                    vy = np.copysign(
                            dy / dt,
                            data['latitude'][idx]
                            - data['latitude'][indexes[0] - 1])
                speed = np.sqrt(dx ** 2 + dy ** 2) / dt
                bearing = (360. + np.rad2deg(np.arctan2(vx, vy))) % 360
                data['speed'].append(int(np.around(speed, 0)))
                data['direction'].append(
                        int(np.around(bearing, 0)))
        return data


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
                    msg = "%s, Retrying in %d seconds..." % (str(e), mdelay)
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
        return urllib.request.urlopen(url)

    res = request_url()
    df = read_csv(
            StringIO("".join([_.decode('utf-8') for _ in res])),
            header=None,
    )
    name = f"{storm_id[:-4].upper():>10}"
    year = f"{storm_id[-4:]:>5}"
    entry = df[(df[0].isin([name]) & df[8].isin([year]))]
    if len(entry) == 0:
        return None
    else:
        return entry[20].tolist()[0].strip()
