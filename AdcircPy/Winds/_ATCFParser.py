import abc
from datetime import datetime
import urllib
import io
import gzip
import numpy as np


class _ATCFParser(metaclass=abc.ABCMeta):
    """
    This class is almost exclusively an ATCF parser.
    """

    def __init__(self, storm_id, start_date=None, end_date=None):
        self._storm_id = storm_id
        self._start_date = start_date
        self._end_date = end_date
        self.__init_vars()
        self.__parse_ATCF()
        self.__cleanup_data()

    def __init_vars(self):
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

    def __parse_ATCF(self):
        response = urllib.request.urlopen(self.url)
        compressed_file = io.BytesIO(response.read())
        self._ATCF = gzip.GzipFile(fileobj=compressed_file)
        for i, line in enumerate(self._ATCF):
            line = line.decode('UTF-8').split(',')
            # filter out lines with no isotach data
            _NEQ = int(line[13].strip(' '))
            _SEQ = int(line[14].strip(' '))
            _SWQ = int(line[15].strip(' '))
            _NWQ = int(line[16].strip(' '))
            _check = np.asarray([_NEQ, _SEQ, _SWQ, _NWQ])
            if not np.all(_check == 0):
                self._basin.append(line[0])
                self._storm_number.append(line[1].strip(' '))
                _datetime = line[2].strip(' ')
                _minutes = line[3].strip(' ')
                if _minutes == '':
                    _minutes = '00'
                _datetime = _datetime+_minutes
                self._datetime.append(datetime.strptime(_datetime,
                                                        '%Y%m%d%H%M'))
                self._record_type.append(line[4].strip(' '))
                if 'N' in line[6]:
                    _lat = float(line[6].strip('N '))*.1
                elif 'S' in line:
                    _lat = float(line[6].strip('S '))*-.1
                self._latitude.append(_lat)
                if 'E' in line[7]:
                    _lon = float(line[7].strip('E '))*.1
                elif 'W' in line[7]:
                    _lon = float(line[7].strip('W '))*-.1
                self._longitude.append(_lon)
                self._max_sustained_wind_speed.append(float(
                                                        line[8].strip(' ')))
                self._central_pressure.append(float(line[9].strip(' ')))
                self._development_level.append(line[10].strip(' '))
                self._isotach.append(int(line[11].strip(' ')))
                self._quadrant.append(line[12].strip(' '))
                self._radius_for_NEQ.append(int(line[13].strip(' ')))
                self._radius_for_SEQ.append(int(line[14].strip(' ')))
                self._radius_for_SWQ.append(int(line[15].strip(' ')))
                self._radius_for_NWQ.append(int(line[16].strip(' ')))
                if len(line) > 18:
                    self._background_pressure.append(int(line[17].strip(' ')))
                    self._radius_of_last_closed_isobar.append(int(
                                                        line[18].strip(' ')))
                    self._radius_of_maximum_winds.append(int(
                                                        line[19].strip(' ')))
                    if len(line) > 23:
                        self._name.append(line[27].strip(' '))
                    else:
                        self._name.append('')

                else:
                    self._background_pressure.append(None)
                    self._radius_of_last_closed_isobar.append(None)
                    self._radius_of_maximum_winds.append(None)
                    self._name.append('')

    def __cleanup_data(self):
        if self.start_date is not None or self.end_date is not None:
            if isinstance(self.start_date, datetime) is False:
                raise Exception("start_date must be a datetime.datetime instance.")
            if isinstance(self.end_date, datetime) is False:
                raise Exception("end_date must be a datetime.datetime instance.")
            if self.end_date <= self.start_date:
                raise Exception("start_date must be previous to end_date.")
            if (self.datetime[0] > self.end_date
                    or self.end_date > self.datetime[-1]):
                raise Exception("end_date provided is out of range with the "
                                + "record\n end_date provided is "
                                + "{}".format(self.end_date.strftime(
                                                            '%Y-%m-%d %H:%M'))
                                + " but the ""last available data point is "
                                + "{}".format(self.datetime[-1].strftime(
                                                            '%Y-%m-%d %H:%M')))
            if (self.datetime[0] > self.start_date
                    or self.start_date > self.datetime[-1]):
                raise Exception("start_date provided is out of range with "
                                + "the record")
        else:
            self._start_date = self.datetime[0]
            self._end_date = self.datetime[-1]
        # get the start and end indexes based on dates provided
        _diff = [np.abs((x - self.start_date).total_seconds())
                 for x in self.datetime]
        si = _diff.index(min(_diff))
        _diff = [np.abs((self.end_date - x).total_seconds())
                 for x in self.datetime]
        ei = _diff.index(min(_diff)) + 1
        # crop the dataset to only include dates provided
        self._basin = self._basin[si:ei]
        self._storm_number = self._storm_number[si:ei]
        self._record_type = self._record_type[si:ei]
        self._latitude = self._latitude[si:ei]
        self._longitude = self._longitude[si:ei]
        self._datetime = self._datetime[si:ei]
        self._max_sustained_wind_speed = self._max_sustained_wind_speed[si:ei]
        self._central_pressure = self._central_pressure[si:ei]
        self._development_level = self._development_level[si:ei]
        self._isotach = self._isotach[si:ei]
        self._quadrant = self._quadrant[si:ei]
        self._radius_for_NEQ = self._radius_for_NEQ[si:ei]
        self._radius_for_SEQ = self._radius_for_SEQ[si:ei]
        self._radius_for_SWQ = self._radius_for_SWQ[si:ei]
        self._radius_for_NWQ = self._radius_for_NWQ[si:ei]
        self._background_pressure = self._background_pressure[si:ei]
        self._radius_of_last_closed_isobar = \
            self._radius_of_last_closed_isobar[si:ei]
        self._radius_of_maximum_winds = self._radius_of_maximum_winds[si:ei]
        self._name = self._name[si:ei]

    @property
    @abc.abstractmethod
    def url(self):
        return self._url

    @property
    def start_date(self):
        return self._start_date

    @property
    def end_date(self):
        return self._end_date

    @property
    def datetime(self):
        return self._datetime

    @property
    def storm_id(self):
        return self._storm_id

    @property
    def longitude(self):
        return np.asarray(self._longitude)

    @property
    def latitude(self):
        return np.asarray(self._latitude)
