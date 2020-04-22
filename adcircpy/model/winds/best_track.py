import numpy as np
import urllib.request
import io
import gzip
from datetime import datetime
from pandas import DataFrame
import pathlib
from shapely.geometry import Point, Polygon
import utm
from haversine import haversine
from pyproj import Proj
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from adcircpy.model.winds.wind_forcing import WindForcing

# import os
# from pathlib import Path
# import zipfile
# from adcircpy.lib._get_cache_directory import _get_cache_directory


class BestTrackForcing(WindForcing):

    def __init__(self, storm_id, start_date=None, end_date=None, dst_crs=None):
        self._storm_id = storm_id
        self._start_date = start_date
        self._end_date = end_date
        self._dst_crs = dst_crs

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
            radii = 1852.*radii  # convert to meters
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
            U = self.speed.iloc[i]*np.sin(np.deg2rad(self.direction.iloc[i]))
            V = self.speed.iloc[i]*np.cos(np.deg2rad(self.direction.iloc[i]))
            ax.quiver(
                self.longitude.iloc[i], self.latitude.iloc[i], U, V, **kwargs)
            ax.annotate(
                self.df['datetime'].iloc[i],
                (self.longitude.iloc[i], self.latitude.iloc[i])
                )
        if show:
            ax.axis('scaled')
            plt.show()

    def write(self, path, overwrite=False):
        path = pathlib.Path(path)
        if path.is_file() and not overwrite:
            raise Exception(
                'Files exist, set overwrite=True to allow overwrite.')
        with open(path, 'w') as f:
            f.write(self.fort22)

    @property
    def storm_id(self):
        self._storm_id

    @property
    def start_date(self):
        return self._start_date

    @property
    def end_date(self):
        return self._end_date

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
    def fort22(self):
        record_number = self._generate_record_numbers()
        fort22 = ''
        for i, (_, row) in enumerate(self.df.iterrows()):
            fort22 += "{:<2},".format(row["basin"])
            fort22 += "{:>3},".format(row["storm_number"])
            fort22 += "{:>11},".format(row["datetime"].strftime('%Y%m%d%H'))
            fort22 += "{:3},".format("")
            fort22 += "{:>5},".format(row["record_type"])
            fort22 += "{:>4},".format(int((row["datetime"]-self.start_date)
                                          .total_seconds()/3600))
            if row["latitude"] >= 0:
                fort22 += "{:>4}N,".format(int(row["latitude"]/.1))
            else:
                fort22 += "{:>4}S,".format(int(row["latitude"]/-.1))
            if row["longitude"] >= 0:
                fort22 += "{:>5}E,".format(int(row["longitude"]/.1))
            else:
                fort22 += "{:>5}W,".format(int(row["longitude"]/-.1))
            fort22 += "{:>4},".format(int(row["max_sustained_wind_speed"]))
            fort22 += "{:>5},".format(int(row["central_pressure"]))
            fort22 += "{:>3},".format(row["development_level"])
            fort22 += "{:>4},".format(int(row["isotach"]))
            fort22 += "{:>4},".format(row["quadrant"])
            fort22 += "{:>5},".format(int(row["radius_for_NEQ"]))
            fort22 += "{:>5},".format(int(row["radius_for_SEQ"]))
            fort22 += "{:>5},".format(int(row["radius_for_SWQ"]))
            fort22 += "{:>5},".format(int(row["radius_for_NWQ"]))
            if row["background_pressure"] is None:
                row["background_pressure"] = \
                    self.df["background_pressure"].iloc[i-1]
            if (row["background_pressure"] <= row["central_pressure"]
                    and 1013 > row["central_pressure"]):
                fort22 += "{:>5},".format(1013)
            elif (row["background_pressure"] <= row["central_pressure"]
                  and 1013 <= row["central_pressure"]):
                fort22 += "{:>5},".format(int(row["central_pressure"]+1))
            else:
                fort22 += "{:>5},".format(int(row["background_pressure"]))
            fort22 += "{:>5},".format(int(
                                        row["radius_of_last_closed_isobar"]))
            fort22 += "{:>4},".format(int(row["radius_of_maximum_winds"]))
            fort22 += "{:>5},".format('')  # gust
            fort22 += "{:>4},".format('')  # eye
            fort22 += "{:>4},".format('')  # subregion
            fort22 += "{:>4},".format('')  # maxseas
            fort22 += "{:>4},".format('')  # initials
            fort22 += "{:>3},".format(row["direction"])
            fort22 += "{:>4},".format(row["speed"])
            fort22 += "{:^12},".format(row["name"])
            # from this point forwards it's all aswip
            fort22 += "{:>4},".format(record_number[i])
            fort22 += "\n"
        return fort22

    @property
    def _storm_id(self):
        storm_id = f"{self.basin}"
        storm_id += f"{self.storm_number}"
        storm_id += f"{self.year}"
        return storm_id

    @property
    def _start_date(self):
        return self.__start_date

    @property
    def _end_date(self):
        return self.__end_date

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
                _datetime = _datetime+_minutes
                data['datetime'].append(
                    datetime.strptime(_datetime, '%Y%m%d%H%M'))
                data['record_type'].append(line[4].strip(' '))
                if 'N' in line[6]:
                    _lat = float(line[6].strip('N '))*.1
                elif 'S' in line:
                    _lat = float(line[6].strip('S '))*-.1
                data['latitude'].append(_lat)
                if 'E' in line[7]:
                    _lon = float(line[7].strip('E '))*.1
                elif 'W' in line[7]:
                    _lon = float(line[7].strip('W '))*-.1
                data['longitude'].append(_lon)
                data['max_sustained_wind_speed'].append(
                    float(line[8].strip(' ')))
                data['central_pressure'].append(
                    float(line[9].strip(' ')))
                data['development_level'].append(
                    line[10].strip(' '))
                data['isotach'].append(int(line[11].strip(' ')))
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

    def _generate_record_numbers(self):
        record_number = [1]
        for i in range(1, len(self.datetime)):
            if self.datetime.iloc[i] == self.datetime.iloc[i-1]:
                record_number.append(record_number[-1])
            else:
                record_number.append(record_number[-1] + 1)
        return record_number

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
                if indexes[-1]+1 < len(data['datetime']):
                    dt = ((data['datetime'][indexes[-1]+1]
                           - data['datetime'][idx])
                          .total_seconds()/(60.*60.))
                    dx = haversine(
                        (data['latitude'][idx],
                         data['longitude'][indexes[-1]+1]),
                        (data['latitude'][idx],
                         data['longitude'][idx]), unit='nmi')
                    dy = haversine(
                        (data['latitude'][indexes[-1]+1],
                         data['longitude'][idx]),
                        (data['latitude'][idx],
                         data['longitude'][idx]), unit='nmi')
                    vx = np.copysign(
                        dx/dt,
                        data['longitude'][indexes[-1]+1]
                        - data['longitude'][idx])
                    vy = np.copysign(
                        dy/dt,
                        data['latitude'][indexes[-1]+1]
                        - data['latitude'][idx])
                else:
                    dt = ((data['datetime'][idx]
                          - data['datetime'][indexes[0]-1])
                          .total_seconds()/(60.*60.))
                    dx = haversine(
                        (data['latitude'][idx],
                         data['longitude'][indexes[0]-1]),
                        (data['latitude'][idx],
                         data['longitude'][idx]), unit='nmi')
                    dy = haversine(
                        (data['latitude'][indexes[0]-1],
                         data['longitude'][idx]),
                        (data['latitude'][idx],
                         data['longitude'][idx]), unit='nmi')
                    vx = np.copysign(
                        dx/dt,
                        data['longitude'][idx]
                        - data['longitude'][indexes[0]-1])
                    vy = np.copysign(
                        dy/dt,
                        data['latitude'][idx]
                        - data['latitude'][indexes[0]-1])
                speed = np.sqrt(dx**2+dy**2)/dt
                bearing = (360. + np.rad2deg(np.arctan2(vx, vy))) % 360
                data['speed'].append(int(np.around(speed, 0)))
                data['direction'].append(
                    int(np.around(bearing, 0)))
        return data

    def transform_to(self, crs):
        pass

    @property
    def NWS(self):
        try:
            return self.__NWS
        except AttributeError:
            return 20

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

    @property
    def geofactor(self):
        try:
            return self.__geofactor
        except AttributeError:
            return 1

    @start_date.setter
    def start_date(self, start_date):
        self._start_date = start_date

    @end_date.setter
    def end_date(self, end_date):
        self._end_date = end_date

    @NWS.setter
    def NWS(self, NWS):
        assert NWS in [19, 20]
        self.__NWS = int(NWS)

    @BLADj.setter
    def BLADj(self, BLADj):
        BLADj = float(BLADj)
        assert BLADj >= 0 and BLADj <= 1
        self.__BLADj = BLADj

    @geofactor.setter
    def geofactor(self, geofactor):
        geofactor = float(geofactor)
        assert geofactor >= 0 and geofactor <= 1
        self.__geofactor = geofactor

    @property
    def _file_end_date(self):
        unique_dates = np.unique(self._df['datetime'])
        for date in unique_dates:
            if date >= self.end_date:
                return date

    @_storm_id.setter
    def _storm_id(self, storm_id):
        url = 'ftp://ftp.nhc.noaa.gov/atcf/archive/'
        url += storm_id[4:]
        url += '/b'
        url += storm_id[0:2].lower()
        url += storm_id[2:]
        url += '.dat.gz'
        try:
            response = urllib.request.urlopen(url)
        except urllib.error.URLError:
            raise NameError(
                f'Did not find storm with id {storm_id}.'
                + f'Submitted URL was {url}.')
        self.__atcf = io.BytesIO(response.read())

    @_start_date.setter
    def _start_date(self, start_date):
        if start_date is not None:
            assert isinstance(start_date, datetime)
        else:
            start_date = self._df['datetime'].iloc[0]
        msg = f"start_date must be >= {self._df['datetime'].iloc[0]} "
        msg += f"and <{self._df['datetime'].iloc[-1]}"
        assert start_date >= self._df['datetime'].iloc[0] \
            and start_date < self._df['datetime'].iloc[-1], msg
        self.__start_date = start_date

    @_end_date.setter
    def _end_date(self, end_date):
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


    # def is_constant_dt(self):
    #     """ returns True if output data is at constant dt """
    #     if all(v.total_seconds() == 0
    #             for v in np.diff(np.unique(self.datetime), n=2)):
    #         return True
    #     else:
    #         return False

    # def remove_TS(self, which='leading'):
    #     assert which in ['leading', 'all']
    #     if which == 'all':
    #         mask = np.logical_or(
    #             self._container['development_level'] == 'TS',
    #             self.__get_active_mask())
    #     else:
    #         mask = self.__get_active_mask()
    #         i = 0
    #         while self._container['development_level'][i] != 'HU':
    #             mask[i] = True
    #             i += 1
    #     self.__set_active_mask(mask)

    # def remove_EX(self, which='trailing'):
    #     assert which in ['trailing', 'all']
    #     if which == 'all':
    #         mask = np.logical_or(
    #             self._container['development_level'] == 'EX',
    #             self.__get_active_mask())
    #     else:
    #         mask = np.flip(self.__get_active_mask())
    #         i = 0
    #         while np.flip(self._container['development_level'])[i] != 'HU':
    #             mask[i] = True
    #             i += 1
    #     self.__set_active_mask(mask)

    # def only_HU(self):
    #     mask = np.full(self._container['datetime'].shape, False)
    #     mask[np.where(self._container['development_level'] == 'HU')] = True
    #     self.__set_active_mask(mask)

    # def remove_non_six_hourly(self):
    #     mask = np.full(self._container['datetime'].shape, False)
    #     for i, datetime in enumerate(self._container['datetime']):
    #         if datetime.hour not in [0, 6, 12, 18]:
    #             mask[i] = True
    #     mask = np.logical_or(mask, self.__get_active_maks())
    #     self.__set_active_mask(mask)

    # def reset_filters(self):
    #     self.__reset_active_mask()

    # def dump(self, path, overwrite=False):
    #     if os.path.isfile(path) and not overwrite:
    #         raise Exception(
    #             'Files exist, set overwrite=True to allow overwrite.')
    #     with open(str(path), 'w') as f:
    #         f.write(self.fort22)

    # def plot_track(self):
    #     coastlines = self.__fetch_coastline()
    #     for coastline in coastlines:
    #         plt.plot(coastline[:, 0], coastline[:, 1], color='k')
    #     for i in range(len(self.speed)):
    #         # when dealing with nautical degrees, U is sine and V is cosine.
    #         U = self.speed[i]*np.sin(np.deg2rad(self.direction[i]))
    #         V = self.speed[i]*np.cos(np.deg2rad(self.direction[i]))
    #         plt.quiver(self.longitude[i], self.latitude[i], U, V)
    #     plt.gca().axis('scaled')
    #     plt.gca().set_xlim([np.min(self.longitude) - 5.,
    #                         np.max(self.longitude) + 5.])
    #     plt.gca().set_ylim([np.min(self.latitude) - 5.,
    #                         np.max(self.latitude) + 5.])
    #     plt.show()

    # def __fetch_coastline(self):
    #     cache = _get_cache_directory()
    #     file_path = str(Path(
    #         cache + '/ne_110m_coastline/ne_110m_coastline.shp'))
    #     if not os.path.isfile(file_path):
    #         zip_file_path = str(Path(cache + '/ne_110m_coastline.zip'))
    #         url = 'http://naciscdn.org/naturalearth/'
    #         url += '110m/physical/ne_110m_coastline.zip'
    #         urllib.request.urlretrieve(url, zip_file_path)
    #         _zip = zipfile.ZipFile(zip_file_path)
    #         outdir = str(Path(cache + '/ne_110m_coastline'))
    #         os.makedirs(outdir, exist_ok=True)
    #         for name in _zip.namelist():
    #             data = _zip.read(name)
    #             outfile = os.path.join(outdir, name)
    #             f = open(outfile, 'wb')
    #             f.write(data)
    #             f.close()
    #         os.remove(zip_file_path)
    #     DataSource = ogr.Open(file_path)
    #     Layer = DataSource.GetLayer()
    #     vertices_collection = list()
    #     for feature in Layer:
    #         Geometry = feature.GetGeometryRef()
    #         points = Geometry.GetPoints()
    #         vertices_collection.append(np.asarray(points))
    #     return vertices_collection

    # def __parse_ATCF(self):
    #     for i, line in enumerate(self.__ATCF):
    #         line = line.decode('UTF-8').split(',')
    #         # filter out lines with no isotach data
    #         _NEQ = int(line[13].strip(' '))
    #         _SEQ = int(line[14].strip(' '))
    #         _SWQ = int(line[15].strip(' '))
    #         _NWQ = int(line[16].strip(' '))
    #         _check = np.asarray([_NEQ, _SEQ, _SWQ, _NWQ])
    #         if not np.all(_check == 0):
    #             self._container['basin'].append(line[0])
    #             self._container['storm_number'].append(line[1].strip(' '))
    #             _datetime = line[2].strip(' ')
    #             _minutes = line[3].strip(' ')
    #             if _minutes == '':
    #                 _minutes = '00'
    #             _datetime = _datetime+_minutes
    #             self._container['datetime'].append(
    #                 datetime.strptime(_datetime, '%Y%m%d%H%M'))
    #             self._container['record_type'].append(line[4].strip(' '))
    #             if 'N' in line[6]:
    #                 _lat = float(line[6].strip('N '))*.1
    #             elif 'S' in line:
    #                 _lat = float(line[6].strip('S '))*-.1
    #             self._container['latitude'].append(_lat)
    #             if 'E' in line[7]:
    #                 _lon = float(line[7].strip('E '))*.1
    #             elif 'W' in line[7]:
    #                 _lon = float(line[7].strip('W '))*-.1
    #             self._container['longitude'].append(_lon)
    #             self._container['max_sustained_wind_speed'].append(
    #                 float(line[8].strip(' ')))
    #             self._container['central_pressure'].append(
    #                 float(line[9].strip(' ')))
    #             self._container['development_level'].append(
    #                 line[10].strip(' '))
    #             self._container['isotach'].append(int(line[11].strip(' ')))
    #             self._container['quadrant'].append(line[12].strip(' '))
    #             self._container['radius_for_NEQ'].append(
    #                 int(line[13].strip(' ')))
    #             self._container['radius_for_SEQ'].append(
    #                 int(line[14].strip(' ')))
    #             self._container['radius_for_SWQ'].append(
    #                 int(line[15].strip(' ')))
    #             self._container['radius_for_NWQ'].append(
    #                 int(line[16].strip(' ')))
    #             if len(line) > 18:
    #                 self._container['background_pressure'].append(
    #                     int(line[17].strip(' ')))
    #                 self._container['radius_of_last_closed_isobar'].append(
    #                     int(line[18].strip(' ')))
    #                 self._container['radius_of_maximum_winds'].append(
    #                     int(line[19].strip(' ')))
    #                 if len(line) > 23:
    #                     self._container['name'].append(line[27].strip(' '))
    #                 else:
    #                     self._container['name'].append('')
    #             else:
    #                 self._container['background_pressure'].append(
    #                     self._container['background_pressure'][-1])
    #                 self._container['radius_of_last_closed_isobar'].append(
    #                     self._container['radius_of_last_closed_isobar'][-1])
    #                 self._container['radius_of_maximum_winds'].append(
    #                     self._container['radius_of_maximum_winds'][-1])
    #                 self._container['name'].append('')
    #     for key, values in self.__container.items():
    #         self._container[key] = np.asarray(values)

    # def __get_active_mask(self):
    #     try:
    #         return self.__active_mask
    #     except AttributeError:
    #         self.__reset_active_mask()
    #         return self.__active_mask

    # def __set_active_mask(self, mask):
    #     self.__active_mask = mask

    # def __reset_active_mask(self):
    #     mask = np.full(len(self._container['datetime']))
    #     self.__active_mask = np.ma.asarray(self._container['datetime']).mask

    # def __init_velocity(self):
    #     self._container['speed'] = list()
    #     self._container['direction'] = list()
    #     zone = utm.from_latlon(
    #         self._container['latitude'][0],
    #         self._container['longitude'][0])[2]
    #     utm_proj = pyproj.Proj(proj='utm', zone=zone)
    #     x, y = utm_proj(
    #         self._container['longitude'], self._container['latitude'])
    #     unique_datetimes = np.unique(self._container['datetime'])
    #     for i, datetime in enumerate(unique_datetimes):
    #         indexes, = np.where(
    #             np.asarray(self._container['datetime']) == datetime)
    #         for idx in indexes:
    #             if indexes[-1]+1 < len(self._container['datetime']):
    #                 dt = ((self._container['datetime'][indexes[-1]+1]
    #                        - self._container['datetime'][idx])
    #                       .total_seconds()/(60.*60.))
    #                 dx = haversine(
    #                     (self._container['latitude'][idx],
    #                      self._container['longitude'][indexes[-1]+1]),
    #                     (self._container['latitude'][idx],
    #                      self._container['longitude'][idx]), unit='nmi')
    #                 dy = haversine(
    #                     (self._container['latitude'][indexes[-1]+1],
    #                      self._container['longitude'][idx]),
    #                     (self._container['latitude'][idx],
    #                      self._container['longitude'][idx]), unit='nmi')
    #                 vx = np.copysign(
    #                     dx/dt,
    #                     self._container['longitude'][indexes[-1]+1]
    #                     - self._container['longitude'][idx])
    #                 vy = np.copysign(
    #                     dy/dt,
    #                     self._container['latitude'][indexes[-1]+1]
    #                     - self._container['latitude'][idx])
    #             else:
    #                 dt = ((self._container['datetime'][idx]
    #                       - self._container['datetime'][indexes[0]-1])
    #                       .total_seconds()/(60.*60.))
    #                 dx = haversine(
    #                     (self._container['latitude'][idx],
    #                      self._container['longitude'][indexes[0]-1]),
    #                     (self._container['latitude'][idx],
    #                      self._container['longitude'][idx]), unit='nmi')
    #                 dy = haversine(
    #                     (self._container['latitude'][indexes[0]-1],
    #                      self._container['longitude'][idx]),
    #                     (self._container['latitude'][idx],
    #                      self._container['longitude'][idx]), unit='nmi')
    #                 vx = np.copysign(
    #                     dx/dt,
    #                     self._container['longitude'][idx]
    #                     - self._container['longitude'][indexes[0]-1])
    #                 vy = np.copysign(
    #                     dy/dt,
    #                     self._container['latitude'][idx]
    #                     - self._container['latitude'][indexes[0]-1])
    #             speed = np.sqrt(dx**2+dy**2)/dt
    #             bearing = (360. + np.rad2deg(np.arctan2(vx, vy))) % 360
    #             self._container['speed'].append(int(np.around(speed, 0)))
    #             self._container['direction'].append(
    #                 int(np.around(bearing, 0)))

    # @property
    # def storm_id(self):
    #     return self.__storm_id



    # @property
    # def url(self):
    #     url = 'ftp://ftp.nhc.noaa.gov/atcf/archive/'
    #     url += self.storm_id[4:]
    #     url += '/b'
    #     url += self.storm_id[0:2].lower()
    #     url += self.storm_id[2:]
    #     url += '.dat.gz'
    #     return url

    # @property
    # def start_date(self):
    #     return self._container['datetime'][self._active_indexes][0]

    # @property
    # def end_date(self):
    #     return self._container['datetime'][self._active_indexes][-1]

    # @property
    # def basin(self):
    #     return self._container['basin'][self._active_indexes]

    # @property
    # def storm_number(self):
    #     return self._container['storm_number'][self._active_indexes]

    # @property
    # def datetime(self):
    #     return self._container['datetime'][self._active_indexes]

    # @property
    # def record_type(self):
    #     return self._container['record_type'][self._active_indexes]

    # @property
    # def latitude(self):
    #     return self._container['latitude'][self._active_indexes]

    # @property
    # def longitude(self):
    #     return self._container['longitude'][self._active_indexes]

    # @property
    # def max_sustained_wind_speed(self):
    #     return self._container[
    #         'max_sustained_wind_speed'][self._active_indexes]

    # @property
    # def central_pressure(self):
    #     return self._container['central_pressure'][self._active_indexes]

    # @property
    # def development_level(self):
    #     return self._container['development_level'][self._active_indexes]

    # @property
    # def isotach(self):
    #     return self._container['isotach'][self._active_indexes]

    # @property
    # def quadrant(self):
    #     return self._container['quadrant'][self._active_indexes]

    # @property
    # def radius_for_NEQ(self):
    #     return self._container['radius_for_NEQ'][self._active_indexes]

    # @property
    # def radius_for_SEQ(self):
    #     return self._container['radius_for_SEQ'][self._active_indexes]

    # @property
    # def radius_for_SWQ(self):
    #     return self._container['radius_for_SWQ'][self._active_indexes]

    # @property
    # def radius_for_NWQ(self):
    #     return self._container['radius_for_NWQ'][self._active_indexes]

    # @property
    # def background_pressure(self):
    #     return self._container['background_pressure'][self._active_indexes]

    # @property
    # def radius_of_last_closed_isobar(self):
    #     return self._container[
    #         'radius_of_last_closed_isobar'][self._active_indexes]

    # @property
    # def radius_of_maximum_winds(self):
    #     return self._container[
    #         'radius_of_maximum_winds'][self._active_indexes]

    # @property
    # def name(self):
    #     return self._container['name'][self._active_indexes]

    # @property
    # def record_number(self):
    #     unique_datetimes = np.unique(self.datetime)
    #     record_number = np.empty(len(self.datetime))
    #     for i, datetime in enumerate(unique_datetimes):
    #         indexes = np.where(np.asarray(self.datetime) == datetime)
    #         for idx in indexes:
    #             record_number[idx] = i+1
    #     return [int(x) for x in record_number]

    # @property
    # def speed(self):
    #     if 'speed' not in self.__container.keys():
    #         self.__init_velocity()
    #     return self._container['speed']

    # @property
    # def direction(self):
    #     if 'direction' not in self.__container.keys():
    #         self.__init_velocity()
    #     return self._container['direction']

    # @property
    # def _container(self):
    #     try:
    #         return self.__container
    #     except AttributeError:
    #         self.__container = dict()
    #         self.__container['basin'] = list()
    #         self.__container['storm_number'] = list()
    #         self.__container['datetime'] = list()
    #         self.__container['record_type'] = list()
    #         self.__container['latitude'] = list()
    #         self.__container['longitude'] = list()
    #         self.__container['max_sustained_wind_speed'] = list()
    #         self.__container['central_pressure'] = list()
    #         self.__container['development_level'] = list()
    #         self.__container['isotach'] = list()
    #         self.__container['quadrant'] = list()
    #         self.__container['radius_for_NEQ'] = list()
    #         self.__container['radius_for_SEQ'] = list()
    #         self.__container['radius_for_SWQ'] = list()
    #         self.__container['radius_for_NWQ'] = list()
    #         self.__container['background_pressure'] = list()
    #         self.__container['radius_of_last_closed_isobar'] = list()
    #         self.__container['radius_of_maximum_winds'] = list()
    #         self.__container['name'] = list()
    #         return self.__container

  
    # @start_date.setter
    # def start_date(self, start_date):
    #     assert isinstance(start_date, datetime)
    #     mask = np.logical_or(
    #         self._container['datetime'] < start_date,
    #         self.__get_active_mask())
    #     self.__set_active_mask(mask)

    # @end_date.setter
    # def end_date(self, end_date):
    #     assert isinstance(end_date, datetime)
    #     mask = np.logical_or(
    #         self._container['datetime'] > end_date,
    #         self.__get_active_mask())
    #     self.__set_active_mask(mask)

    # @property
    # def _active_indexes(self):
    #     return np.where(~self.__get_active_mask())[0]

    # @storm_id.setter
    # def storm_id(self, storm_id):
    #     self.__storm_id = storm_id
    #     self.__fetch_ATCF()
    #     self.__parse_ATCF()


