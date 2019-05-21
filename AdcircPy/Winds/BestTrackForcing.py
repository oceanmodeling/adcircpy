import urllib
import io
import gzip
import numpy as np
from datetime import datetime
import utm
from haversine import haversine
import pyproj
from pathlib import Path

# local imports
from AdcircPy.Winds import _WindForcing

# unittest imports
import unittest


class BestTrackForcing(_WindForcing):
    def __init__(self, storm_id, start_date=None, end_date=None):
        self._storm_id = storm_id
        self._basin = list()
        self._storm_number = list()
        self._datetime = list()
        self._record_type = list()
        self._latitude = list()
        self._longitude = list()
        self._max_sustained_wind_speed = list()
        self._central_pressure = list()
        self._development_level = list()
        self._isotach = list()
        self._quadrant = list()
        self._radius_for_NEQ = list()
        self._radius_for_SEQ = list()
        self._radius_for_SWQ = list()
        self._radius_for_NWQ = list()
        self._background_pressure = list()
        self._radius_of_last_closed_isobar = list()
        self._radius_of_maximum_winds = list()
        self._name = list()
        self._speed = list()
        self._direction = list()
        self.__set_url()
        self.__set_ATCF()
        self.__parse_ATCF()
        self._start_date = start_date
        self._end_date = end_date
        self.__cleanup_data()
        self.__init_speed_and_direction()
        self.__set_record_number()

    def dump(self, path, filename='fort.22.best_track'):
        path = Path(path + '/' + filename)
        fort22 = self.get_fort22()
        with open(str(path), 'w') as f:
            f.write(fort22)

    def get_fort22(self):
        fort22 = ''
        for i in range(len(self.datetime)):
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

    def __set_url(self):
        url = 'ftp://ftp.nhc.noaa.gov/atcf/archive/' \
              + self.storm_id[4:] + '/b' + self.storm_id[0:2].lower() \
              + self.storm_id[2:] + '.dat.gz'
        self.__url = url

    def __set_ATCF(self):
        try:
            response = urllib.request.urlopen(self.url)
        except urllib.error.URLError:
            raise NameError(
                'Did not find storm with id '
                + '{}. '.format(self.storm_id)
                + 'Submitted URL was {}'.format(self.url))

        compressed_file = io.BytesIO(response.read())
        self.__ATCF = gzip.GzipFile(fileobj=compressed_file)

    def __parse_ATCF(self):
        for i, line in enumerate(self. ATCF):
            line = line.decode('UTF-8').split(',')
            # filter out lines with no isotach data
            _NEQ = int(line[13].strip(' '))
            _SEQ = int(line[14].strip(' '))
            _SWQ = int(line[15].strip(' '))
            _NWQ = int(line[16].strip(' '))
            _check = np.asarray([_NEQ, _SEQ, _SWQ, _NWQ])
            if not np.all(_check == 0):
                self.basin.append(line[0])
                self.storm_number.append(line[1].strip(' '))
                _datetime = line[2].strip(' ')
                _minutes = line[3].strip(' ')
                if _minutes == '':
                    _minutes = '00'
                _datetime = _datetime+_minutes
                self.datetime.append(
                    datetime.strptime(_datetime, '%Y%m%d%H%M'))
                self.record_type.append(line[4].strip(' '))
                if 'N' in line[6]:
                    _lat = float(line[6].strip('N '))*.1
                elif 'S' in line:
                    _lat = float(line[6].strip('S '))*-.1
                self.latitude.append(_lat)
                if 'E' in line[7]:
                    _lon = float(line[7].strip('E '))*.1
                elif 'W' in line[7]:
                    _lon = float(line[7].strip('W '))*-.1
                self.longitude.append(_lon)
                self.max_sustained_wind_speed.append(float(line[8].strip(' ')))
                self.central_pressure.append(float(line[9].strip(' ')))
                self.development_level.append(line[10].strip(' '))
                self.isotach.append(int(line[11].strip(' ')))
                self.quadrant.append(line[12].strip(' '))
                self.radius_for_NEQ.append(int(line[13].strip(' ')))
                self.radius_for_SEQ.append(int(line[14].strip(' ')))
                self.radius_for_SWQ.append(int(line[15].strip(' ')))
                self.radius_for_NWQ.append(int(line[16].strip(' ')))
                if len(line) > 18:
                    self.background_pressure.append(int(line[17].strip(' ')))
                    self.radius_of_last_closed_isobar.append(int(
                                                        line[18].strip(' ')))
                    self.radius_of_maximum_winds.append(int(
                                                        line[19].strip(' ')))
                    if len(line) > 23:
                        self.name.append(line[27].strip(' '))
                    else:
                        self.name.append('')
                else:
                    self.background_pressure.append(
                                                self.background_pressure[-1])
                    self.radius_of_last_closed_isobar.append(
                        self.radius_of_last_closed_isobar[-1])
                    self.radius_of_maximum_winds.append(
                        self.radius_of_maximum_winds[-1])
                    self.name.append('')

    def __cleanup_data(self):
        # get the start and end indexes based on dates provided
        _diff = [np.abs((x - self.start_date).total_seconds())
                 for x in self.datetime]
        si = _diff.index(min(_diff))
        _diff = [np.abs((self.end_date - x).total_seconds())
                 for x in self.datetime]
        ei = _diff.index(min(_diff)) + 1
        # crop the dataset to only include dates provided
        self.__basin = self.basin[si:ei]
        self.__storm_number = self.storm_number[si:ei]
        self.__record_type = self.record_type[si:ei]
        self.__latitude = self.latitude[si:ei]
        self.__longitude = self.longitude[si:ei]
        self.__datetime = self.datetime[si:ei]
        self.__max_sustained_wind_speed = self.max_sustained_wind_speed[si:ei]
        self.__central_pressure = self.central_pressure[si:ei]
        self.__development_level = self.development_level[si:ei]
        self.__isotach = self.isotach[si:ei]
        self.__quadrant = self.quadrant[si:ei]
        self.__radius_for_NEQ = self.radius_for_NEQ[si:ei]
        self.__radius_for_SEQ = self.radius_for_SEQ[si:ei]
        self.__radius_for_SWQ = self.radius_for_SWQ[si:ei]
        self.__radius_for_NWQ = self.radius_for_NWQ[si:ei]
        self.__background_pressure = self.background_pressure[si:ei]
        self.__radius_of_last_closed_isobar \
            = self.radius_of_last_closed_isobar[si:ei]
        self.__radius_of_maximum_winds = self.radius_of_maximum_winds[si:ei]
        self.__name = self.name[si:ei]

    def __init_speed_and_direction(self):
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
                self.speed.append(int(np.around(speed, 0)))
                self.direction.append(int(np.around(bearing, 0)))

    def __set_record_number(self):
        unique_datetimes = np.unique(self.datetime)
        record_number = np.empty(len(self.datetime))
        for i, datetime in enumerate(unique_datetimes):
            indexes = np.where(np.asarray(self.datetime) == datetime)
            for idx in indexes:
                record_number[idx] = i+1
        self.__record_number = [int(x) for x in record_number]

    @property
    def storm_id(self):
        return self._storm_id

    @property
    def url(self):
        return self.__url

    @property
    def ATCF(self):
        return self.__ATCF

    @property
    def basin(self):
        return self._basin

    @property
    def storm_number(self):
        return self._storm_number

    @property
    def datetime(self):
        return self._datetime

    @property
    def record_type(self):
        return self._record_type

    @property
    def record_number(self):
        return self.__record_number

    @property
    def latitude(self):
        return self._latitude

    @property
    def longitude(self):
        return self._longitude

    @property
    def max_sustained_wind_speed(self):
        return self._max_sustained_wind_speed

    @property
    def central_pressure(self):
        return self._central_pressure

    @property
    def development_level(self):
        return self._development_level

    @property
    def isotach(self):
        return self._isotach

    @property
    def quadrant(self):
        return self._quadrant

    @property
    def radius_for_NEQ(self):
        return self._radius_for_NEQ

    @property
    def radius_for_SEQ(self):
        return self._radius_for_SEQ

    @property
    def radius_for_SWQ(self):
        return self._radius_for_SWQ

    @property
    def radius_for_NWQ(self):
        return self._radius_for_NWQ

    @property
    def background_pressure(self):
        return self._background_pressure

    @property
    def radius_of_last_closed_isobar(self):
        return self._radius_of_last_closed_isobar

    @property
    def radius_of_maximum_winds(self):
        return self._radius_of_maximum_winds

    @property
    def name(self):
        return self._name

    @property
    def start_date(self):
        return self._start_date

    @property
    def end_date(self):
        return self._end_date

    @property
    def speed(self):
        return self._speed

    @property
    def direction(self):
        return self._direction

    @property
    def _storm_id(self):
        return self.__storm_id

    @property
    def _datetime(self):
        return self.__datetime

    @property
    def _start_date(self):
        return self.__start_date

    @property
    def _end_date(self):
        return self.__end_date

    @property
    def _basin(self):
        return self.__basin

    @property
    def _storm_number(self):
        return self.__storm_number

    @property
    def _datetime(self):
        return self.__datetime

    @property
    def _record_type(self):
        return self.__record_type

    @property
    def _latitude(self):
        return self.__latitude

    @property
    def _longitude(self):
        return self.__longitude

    @property
    def _max_sustained_wind_speed(self):
        return self.__max_sustained_wind_speed

    @property
    def _central_pressure(self):
        return self.__central_pressure

    @property
    def _development_level(self):
        return self.__development_level

    @property
    def _isotach(self):
        return self.__isotach

    @property
    def _quadrant(self):
        return self.__quadrant

    @property
    def _radius_for_NEQ(self):
        return self.__radius_for_NEQ

    @property
    def _radius_for_SEQ(self):
        return self.__radius_for_SEQ

    @property
    def _radius_for_SWQ(self):
        return self.__radius_for_SWQ

    @property
    def _radius_for_NWQ(self):
        return self.__radius_for_NWQ

    @property
    def _background_pressure(self):
        return self.__background_pressure

    @property
    def _radius_of_last_closed_isobar(self):
        return self.__radius_of_last_closed_isobar

    @property
    def _radius_of_maximum_winds(self):
        return self.__radius_of_maximum_winds

    @property
    def _name(self):
        return self.__name

    @property
    def _speed(self):
        return self.__speed

    @property
    def _direction(self):
        return self.__direction

    @_storm_id.setter
    def _storm_id(self, storm_id):
        self.__storm_id = str(storm_id)

    @_basin.setter
    def _basin(self, basin):
        self.__basin = list()

    @_storm_number.setter
    def _storm_number(self, storm_number):
        self.__storm_number = list()

    @_datetime.setter
    def _datetime(self, datetime):
        self.__datetime = list()

    @_record_type.setter
    def _record_type(self, record_type):
        self.__record_type = list()

    @_latitude.setter
    def _latitude(self, latitude):
        self.__latitude = list()

    @_longitude.setter
    def _longitude(self, longitude):
        self.__longitude = list()

    @_max_sustained_wind_speed.setter
    def _max_sustained_wind_speed(self, max_sustained_wind_speed):
        self.__max_sustained_wind_speed = list()

    @_central_pressure.setter
    def _central_pressure(self, central_pressure):
        self.__central_pressure = list()

    @_development_level.setter
    def _development_level(self, development_level):
        self.__development_level = list()

    @_isotach.setter
    def _isotach(self, isotach):
        self.__isotach = list()

    @_quadrant.setter
    def _quadrant(self, quadrant):
        self.__quadrant = list()

    @_radius_for_NEQ.setter
    def _radius_for_NEQ(self, radius_for_NEQ):
        self.__radius_for_NEQ = list()

    @_radius_for_SEQ.setter
    def _radius_for_SEQ(self, radius_for_SEQ):
        self.__radius_for_SEQ = list()

    @_radius_for_SWQ.setter
    def _radius_for_SWQ(self, radius_for_SWQ):
        self.__radius_for_SWQ = list()

    @_radius_for_NWQ.setter
    def _radius_for_NWQ(self, radius_for_NWQ):
        self.__radius_for_NWQ = list()

    @_background_pressure.setter
    def _background_pressure(self, background_pressure):
        self.__background_pressure = list()

    @_radius_of_last_closed_isobar.setter
    def _radius_of_last_closed_isobar(self, radius_of_last_closed_isobar):
        self.__radius_of_last_closed_isobar = list()

    @_radius_of_maximum_winds.setter
    def _radius_of_maximum_winds(self, radius_of_maximum_winds):
        self.__radius_of_maximum_winds = list()

    @_name.setter
    def _name(self, name):
        self.__name = list()

    @_speed.setter
    def _speed(self, speed):
        self.__speed = list()

    @_direction.setter
    def _direction(self, direction):
        self.__direction = list()

    @_start_date.setter
    def _start_date(self, start_date):
        if start_date is not None:
            assert isinstance(start_date, datetime)
            if (self.datetime[0] > start_date
                    or self.start_date > self.datetime[-1]):
                raise Exception("start_date provided is out of range with "
                                + "the record")
        else:
            start_date = self.datetime[0]
        self.__start_date = start_date

    @_end_date.setter
    def _end_date(self, end_date):
        if end_date is not None:
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
        else:
            end_date = self.datetime[-1]
        self.__end_date = end_date


class BestTrackForcingTestCase(unittest.TestCase):

    def test_Sandy(self):
        BTF = BestTrackForcing('AL182012')
        fort22 = BTF.get_fort22()
        print(fort22)

    def test_badname(self):

        BestTrackForcing('badname')
