import urllib.request
import io
import gzip
import os
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import utm
from haversine import haversine
import pyproj
from pathlib import Path
from osgeo import ogr
import zipfile
from adcircpy.lib._WindForcing import _WindForcing
from adcircpy.lib._get_cache_directory import _get_cache_directory


class BestTrackForcing(_WindForcing):

    def __init__(self):
        self.__container = dict()
        self.__container['basin'] = list()
        self.__container['storm_number'] = list()
        self.__container['datetime'] = list()
        self.__container['record_type'] = list()
        self.__container['latitude'] = list()
        self.__container['longitude'] = list()
        self.__container['max_sustained_wind_speed'] = list()
        self.__container['central_pressure'] = list()
        self.__container['development_level'] = list()
        self.__container['isotach'] = list()
        self.__container['quadrant'] = list()
        self.__container['radius_for_NEQ'] = list()
        self.__container['radius_for_SEQ'] = list()
        self.__container['radius_for_SWQ'] = list()
        self.__container['radius_for_NWQ'] = list()
        self.__container['background_pressure'] = list()
        self.__container['radius_of_last_closed_isobar'] = list()
        self.__container['radius_of_maximum_winds'] = list()
        self.__container['name'] = list()

    def clip_data_by_datetime(self, start_date, end_date):
        idx = list(np.where(np.asarray(self.datetime) < start_date)[0])
        for i in reversed(idx):
            for key in list(self.__container.keys()):
                self.__container[key].pop(i)

        idx = list(np.where(np.asarray(self.datetime) > end_date)[0])
        for i in reversed(idx):
            for key in list(self.__container.keys()):
                self.__container[key].pop(i)

    def is_constant_dt(self):
        """ returns True if output data is at constant dt """
        if all(v.total_seconds() == 0
                for v in np.diff(np.unique(self.datetime), n=2)):
            return True
        else:
            return False

    def only_HU(self):
        idx = []
        for i, datetime in enumerate(self.datetime):
            if self.development_level[i] != 'HU':
                idx.append(i)
        for i in reversed(idx):
            for key in list(self.__container.keys()):
                self.__container[key].pop(i)

    def remove_TS(self):
        idx = []
        for i, datetime in enumerate(self.datetime):
            if self.development_level[i] == 'TS':
                idx.append(i)
        for i in reversed(idx):
            for key in list(self.__container.keys()):
                self.__container[key].pop(i)

    def remove_EX(self):
        idx = []
        for i, datetime in enumerate(self.datetime):
            if self.development_level[i] == 'EX':
                idx.append(i)
        for i in reversed(idx):
            for key in list(self.__container.keys()):
                self.__container[key].pop(i)

    def remove_non_six_hourly(self):
        idx = []
        for i, datetime in enumerate(self.datetime):
            if datetime.hour not in [0, 6, 12, 18]:
                idx.append(i)
        for i in reversed(idx):
            for key in list(self.__container.keys()):
                self.__container[key].pop(i)

    def dump(self, path, overwrite=False):
        if os.path.isfile(path) and not overwrite:
            raise Exception(
                'Files exist, set overwrite=True to allow overwrite.')
        with open(str(path), 'w') as f:
            f.write(self.fort22)

    def plot_track(self):
        coastlines = self.__fetch_coastline()
        for coastline in coastlines:
            plt.plot(coastline[:, 0], coastline[:, 1], color='k')
        for i in range(len(self.speed)):
            # when dealing with nautical degrees, U is sine and V is cosine.
            U = self.speed[i]*np.sin(np.deg2rad(self.direction[i]))
            V = self.speed[i]*np.cos(np.deg2rad(self.direction[i]))
            plt.quiver(self.longitude[i], self.latitude[i], U, V)
        plt.gca().axis('scaled')
        plt.gca().set_xlim([np.min(self.longitude) - 5.,
                            np.max(self.longitude) + 5.])
        plt.gca().set_ylim([np.min(self.latitude) - 5.,
                            np.max(self.latitude) + 5.])
        plt.show()

    def __fetch_ATCF(self):
        try:
            response = urllib.request.urlopen(self.url)
        except urllib.error.URLError:
            raise NameError(
                'Did not find storm with id '
                + '{}. '.format(self.storm_id)
                + 'Submitted URL was {}'.format(self.url))
        compressed_file = io.BytesIO(response.read())
        self.__ATCF = gzip.GzipFile(fileobj=compressed_file)

    def __fetch_coastline(self):
        cache = _get_cache_directory()
        file_path = str(Path(
            cache + '/ne_110m_coastline/ne_110m_coastline.shp'))
        if not os.path.isfile(file_path):
            zip_file_path = str(Path(cache + '/ne_110m_coastline.zip'))
            url = 'http://naciscdn.org/naturalearth/'
            url += '110m/physical/ne_110m_coastline.zip'
            urllib.request.urlretrieve(url, zip_file_path)
            _zip = zipfile.ZipFile(zip_file_path)
            outdir = str(Path(cache + '/ne_110m_coastline'))
            os.makedirs(outdir, exist_ok=True)
            for name in _zip.namelist():
                data = _zip.read(name)
                outfile = os.path.join(outdir, name)
                f = open(outfile, 'wb')
                f.write(data)
                f.close()
            os.remove(zip_file_path)
        DataSource = ogr.Open(file_path)
        Layer = DataSource.GetLayer()
        vertices_collection = list()
        for feature in Layer:
            Geometry = feature.GetGeometryRef()
            points = Geometry.GetPoints()
            vertices_collection.append(np.asarray(points))
        return vertices_collection

    def __parse_ATCF(self):
        for i, line in enumerate(self.__ATCF):
            line = line.decode('UTF-8').split(',')
            # filter out lines with no isotach data
            _NEQ = int(line[13].strip(' '))
            _SEQ = int(line[14].strip(' '))
            _SWQ = int(line[15].strip(' '))
            _NWQ = int(line[16].strip(' '))
            _check = np.asarray([_NEQ, _SEQ, _SWQ, _NWQ])
            if not np.all(_check == 0):
                self.__container['basin'].append(line[0])
                self.__container['storm_number'].append(line[1].strip(' '))
                _datetime = line[2].strip(' ')
                _minutes = line[3].strip(' ')
                if _minutes == '':
                    _minutes = '00'
                _datetime = _datetime+_minutes
                self.__container['datetime'].append(
                    datetime.strptime(_datetime, '%Y%m%d%H%M'))
                self.__container['record_type'].append(line[4].strip(' '))
                if 'N' in line[6]:
                    _lat = float(line[6].strip('N '))*.1
                elif 'S' in line:
                    _lat = float(line[6].strip('S '))*-.1
                self.__container['latitude'].append(_lat)
                if 'E' in line[7]:
                    _lon = float(line[7].strip('E '))*.1
                elif 'W' in line[7]:
                    _lon = float(line[7].strip('W '))*-.1
                self.__container['longitude'].append(_lon)
                self.__container['max_sustained_wind_speed'].append(
                    float(line[8].strip(' ')))
                self.__container['central_pressure'].append(
                    float(line[9].strip(' ')))
                self.__container['development_level'].append(
                    line[10].strip(' '))
                self.__container['isotach'].append(int(line[11].strip(' ')))
                self.__container['quadrant'].append(line[12].strip(' '))
                self.__container['radius_for_NEQ'].append(
                    int(line[13].strip(' ')))
                self.__container['radius_for_SEQ'].append(
                    int(line[14].strip(' ')))
                self.__container['radius_for_SWQ'].append(
                    int(line[15].strip(' ')))
                self.__container['radius_for_NWQ'].append(
                    int(line[16].strip(' ')))
                if len(line) > 18:
                    self.__container['background_pressure'].append(
                        int(line[17].strip(' ')))
                    self.__container['radius_of_last_closed_isobar'].append(
                        int(line[18].strip(' ')))
                    self.__container['radius_of_maximum_winds'].append(
                        int(line[19].strip(' ')))
                    if len(line) > 23:
                        self.__container['name'].append(line[27].strip(' '))
                    else:
                        self.__container['name'].append('')
                else:
                    self.__container['background_pressure'].append(
                        self.__container['background_pressure'][-1])
                    self.__container['radius_of_last_closed_isobar'].append(
                        self.__container['radius_of_last_closed_isobar'][-1])
                    self.__container['radius_of_maximum_winds'].append(
                        self.__container['radius_of_maximum_winds'][-1])
                    self.__container['name'].append('')
        self.__init_velocity()

    def __init_velocity(self):
        self.__container['speed'] = list()
        self.__container['direction'] = list()
        zone = utm.from_latlon(self.latitude[0], self.longitude[0])[2]
        utm_proj = pyproj.Proj(proj='utm', zone=zone)
        x, y = utm_proj(self.longitude, self.latitude)
        unique_datetimes = np.unique(self.datetime)
        for i, datetime in enumerate(unique_datetimes):
            indexes, = np.where(np.asarray(self.datetime) == datetime)
            for idx in indexes:
                if indexes[-1]+1 < len(self.datetime):
                    dt = ((self.datetime[indexes[-1]+1] - self.datetime[idx])
                          .total_seconds()/(60.*60.))
                    dx = haversine((self.latitude[idx],
                                    self.longitude[indexes[-1]+1]),
                                   (self.latitude[idx],
                                    self.longitude[idx]), unit='nmi')
                    dy = haversine((self.latitude[indexes[-1]+1],
                                    self.longitude[idx]),
                                   (self.latitude[idx],
                                    self.longitude[idx]), unit='nmi')
                    vx = np.copysign(dx/dt, self.longitude[indexes[-1]+1] -
                                     self.longitude[idx])
                    vy = np.copysign(dy/dt, self.latitude[indexes[-1]+1] -
                                     self.latitude[idx])
                else:
                    dt = ((self.datetime[idx]-self.datetime[indexes[0]-1])
                          .total_seconds()/(60.*60.))
                    dx = haversine((self.latitude[idx],
                                    self.longitude[indexes[0]-1]),
                                   (self.latitude[idx],
                                    self.longitude[idx]), unit='nmi')
                    dy = haversine((self.latitude[indexes[0]-1],
                                    self.longitude[idx]),
                                   (self.latitude[idx],
                                    self.longitude[idx]), unit='nmi')
                    vx = np.copysign(dx/dt, self.longitude[idx] -
                                     self.longitude[indexes[0]-1])
                    vy = np.copysign(dy/dt, self.latitude[idx] -
                                     self.latitude[indexes[0]-1])
                speed = np.sqrt(dx**2+dy**2)/dt
                bearing = (360. + np.rad2deg(np.arctan2(vx, vy))) % 360
                self.__container['speed'].append(int(np.around(speed, 0)))
                self.__container['direction'].append(
                    int(np.around(bearing, 0)))

    @property
    def storm_id(self):
        try:
            return self.__storm_id
        except AttributeError:
            raise AttributeError(
                'Must set storm_id before operation can take place.')

    @property
    def fort22(self):
        fort22 = ''
        for i in self._index_view:
            fort22 += "{:<2},".format(self.basin[i])
            fort22 += "{:>3},".format(self.storm_number[i])
            fort22 += "{:>11},".format(self.datetime[i].strftime('%Y%m%d%H'))
            fort22 += "{:3},".format("")
            fort22 += "{:>5},".format(self.record_type[i])
            fort22 += "{:>4},".format(int((self.datetime[i]-self.start_date)
                                          .total_seconds()/3600))
            if self.latitude[i] >= 0:
                fort22 += "{:>4}N,".format(int(self.latitude[i]/.1))
            else:
                fort22 += "{:>4}S,".format(int(self.latitude[i]/-.1))
            if self.longitude[i] >= 0:
                fort22 += "{:>5}E,".format(int(self.longitude[i]/.1))
            else:
                fort22 += "{:>5}W,".format(int(self.longitude[i]/-.1))
            fort22 += "{:>4},".format(int(self.max_sustained_wind_speed[i]))
            fort22 += "{:>5},".format(int(self.central_pressure[i]))
            fort22 += "{:>3},".format(self.development_level[i])
            fort22 += "{:>4},".format(int(self.isotach[i]))
            fort22 += "{:>4},".format(self.quadrant[i])
            fort22 += "{:>5},".format(int(self.radius_for_NEQ[i]))
            fort22 += "{:>5},".format(int(self.radius_for_SEQ[i]))
            fort22 += "{:>5},".format(int(self.radius_for_SWQ[i]))
            fort22 += "{:>5},".format(int(self.radius_for_NWQ[i]))
            if self.background_pressure is None:
                self.background_pressure[i] = self.background_pressure[i-1]
            if (self.background_pressure[i] <= self.central_pressure[i]
                    and 1013 > self.central_pressure[i]):
                fort22 += "{:>5},".format(1013)
            elif (self.background_pressure[i] <= self.central_pressure[i]
                  and 1013 <= self.central_pressure[i]):
                fort22 += "{:>5},".format(int(self.central_pressure[i]+1))
            else:
                fort22 += "{:>5},".format(int(self.background_pressure[i]))
            fort22 += "{:>5},".format(int(
                                        self.radius_of_last_closed_isobar[i]))
            fort22 += "{:>4},".format(int(self.radius_of_maximum_winds[i]))
            fort22 += "{:>5},".format('')  # gust
            fort22 += "{:>4},".format('')  # eye
            fort22 += "{:>4},".format('')  # subregion
            fort22 += "{:>4},".format('')  # maxseas
            fort22 += "{:>4},".format('')  # initials
            fort22 += "{:>3},".format(self.direction[i])
            fort22 += "{:>4},".format(self.speed[i])
            fort22 += "{:^12},".format(self.name[i])
            # from this point forwards it's all aswip
            fort22 += "{:>4},".format(self.record_number[i])
            fort22 += "\n"
        return fort22

    @property
    def url(self):
        if not hasattr(self, "__url"):
            url = 'ftp://ftp.nhc.noaa.gov/atcf/archive/'
            url += self.storm_id[4:]
            url += '/b'
            url += self.storm_id[0:2].lower()
            url += self.storm_id[2:]
            url += '.dat.gz'
            self.__url = url
        return self.__url

    @property
    def basin(self):
        return self.__container['basin']

    @property
    def storm_number(self):
        return self.__container['storm_number']

    @property
    def datetime(self):
        return self.__container['datetime']

    @property
    def record_type(self):
        return self.__container['record_type']

    @property
    def latitude(self):
        return self.__container['latitude']

    @property
    def longitude(self):
        return self.__container['longitude']

    @property
    def max_sustained_wind_speed(self):
        return self.__container['max_sustained_wind_speed']

    @property
    def central_pressure(self):
        return self.__container['central_pressure']

    @property
    def development_level(self):
        return self.__container['development_level']

    @property
    def isotach(self):
        return self.__container['isotach']

    @property
    def quadrant(self):
        return self.__container['quadrant']

    @property
    def radius_for_NEQ(self):
        return self.__container['radius_for_NEQ']

    @property
    def radius_for_SEQ(self):
        return self.__container['radius_for_SEQ']

    @property
    def radius_for_SWQ(self):
        return self.__container['radius_for_SWQ']

    @property
    def radius_for_NWQ(self):
        return self.__container['radius_for_NWQ']

    @property
    def background_pressure(self):
        return self.__container['background_pressure']

    @property
    def radius_of_last_closed_isobar(self):
        return self.__container['radius_of_last_closed_isobar']

    @property
    def radius_of_maximum_winds(self):
        return self.__container['radius_of_maximum_winds']

    @property
    def name(self):
        return self.__container['name']

    @property
    def record_number(self):
        unique_datetimes = np.unique(self.datetime)
        record_number = np.empty(len(self.datetime))
        for i, datetime in enumerate(unique_datetimes):
            indexes = np.where(np.asarray(self.datetime) == datetime)
            for idx in indexes:
                record_number[idx] = i+1
        return [int(x) for x in record_number]

    @property
    def start_date(self):
        try:
            return self.__start_date
        except AttributeError:
            return self.__container['datetime'][0]

    @property
    def end_date(self):
        try:
            return self.__end_date
        except AttributeError:
            return self.__container['datetime'][-1]

    @property
    def speed(self):
        if 'speed' not in self.__container.keys():
            self.__init_velocity()
        return self.__container['speed']

    @property
    def direction(self):
        if 'direction' not in self.__container.keys():
            self.__init_velocity()
        return self.__container['direction']

    @start_date.setter
    def start_date(self, start_date):
        if isinstance(start_date, str):
            start_date = datetime.strptime(start_date, '%Y%m%d%H')
        assert isinstance(start_date, datetime)
        if (self.datetime[0] > start_date
                or self.start_date > self.datetime[-1]):
            raise Exception("start_date provided is out of range with "
                            + "the record")
        self.clip_data_by_datetime(start_date, self.end_date)

    @end_date.setter
    def end_date(self, end_date):
        self.__parse_ATCF()
        if isinstance(end_date, str):
            end_date = datetime.strptime(end_date, '%Y%m%d%H')
        assert isinstance(end_date, datetime)
        if end_date <= self.start_date:
            raise Exception("start_date must be previous to end_date.")
        if (self.datetime[0] > end_date
                or end_date > self.datetime[-1]):
            raise Exception("end_date provided is out of range with the "
                            + "record\n end_date provided is "
                            + "{}".format(self.end_date.strftime(
                                                        '%Y-%m-%d %H:%M'))
                            + " but the ""last available data point is "
                            + "{}".format(self.datetime[-1].strftime(
                                                        '%Y-%m-%d %H:%M')))
        self.clip_data_by_datetime(self.start_date, end_date)

    @storm_id.setter
    def storm_id(self, storm_id):
        self.__storm_id = storm_id
        self.__fetch_ATCF()
        self.__parse_ATCF()

    @property
    def _index_view(self):
        if not hasattr(self, "__index_view"):
            self._index_view = np.arange(len(self.__container['datetime']))
        return self.__index_view

    @_index_view.setter
    def _index_view(self, indexes):
        self.__index_view = indexes
