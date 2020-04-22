from collections.abc import Mapping
import numpy as np
from datetime import datetime, timedelta
import calendar
import json
import requests


class TidalStations(Mapping):

    def __init__(self):
        self.__storage = dict()

    def __getitem__(self, key):
        return self.__storage[key]

    def __iter__(self):
        return iter(self.__storage)

    def __len__(self):
        return len(self.__storage.keys())

    def add_station(self, station_id, start_date, end_date):
        self.__storage[station_id] = self.__fetch_station_data(
            station_id, start_date, end_date)

    def __fetch_station_data(self, station_id, start_date, end_date):
        responses = list()
        for _start_date, _end_date in self.__get_datetime_segments(
                start_date, end_date):
            params = self.__get_params(station_id, _start_date, _end_date)
            try:
                r = requests.get(self.url, params=params, timeout=10.)
                r.raise_for_status()
            except requests.exceptions.HTTPError as errh:
                print("Http Error:", errh)
            except requests.exceptions.ConnectionError as errc:
                print("Error Connecting:", errc)
            except requests.exceptions.Timeout as errt:
                print("Timeout Error:", errt)
            except requests.exceptions.RequestException as err:
                print("Unknown error.", err)
            responses.append(r)
        data = dict()
        data['datetime'] = list()
        data['values'] = list()
        for i, response in enumerate(responses):
            json_data = json.loads(response.text)
            if 'error' in json_data.keys():
                _start_date, _end_date = list(
                    self.__get_datetime_segments(start_date, end_date))[i]
                data['datetime'].append(_start_date)
                data['values'].append(np.nan)
                data['datetime'].append(_end_date)
                data['values'].append(np.nan)
                continue
            if 'x' not in data.keys():
                data['x'] = float(json_data['metadata']['lon'])
            if 'y' not in data.keys():
                data['y'] = float(json_data['metadata']['lat'])
            if 'name' not in data.keys():
                data['name'] = json_data['metadata']['name']
            for _data in json_data['data']:
                data['datetime'].append(
                    datetime.strptime(_data['t'], '%Y-%m-%d %H:%M'))
                try:
                    data['values'].append(float(_data['v']))
                except ValueError:
                    data['values'].append(np.nan)
        if 'name' not in data.keys():
            data['name'] = ''
        return data

    def __get_params(self, station_id, start_date, end_date):
        params = {}
        params['station'] = station_id
        params['begin_date'] = start_date.strftime('%Y%m%d %H:%M')
        params['end_date'] = end_date.strftime('%Y%m%d %H:%M')
        params['product'] = 'water_level'
        params['datum'] = self.datum
        params['units'] = self.units
        params['time_zone'] = self.time_zone
        params['format'] = 'json'
        params['application'] = 'noaa/nos/csdl/adcircpy'
        return params

    def __get_datetime_segments(self, start_date, end_date):
        """
        https://www.ianwootten.co.uk/2014/07/01/splitting-a-date-range-in-python/
        """

        segments = [(start_date, end_date)]
        interval = 2
        while np.any([(_end_date-_start_date).total_seconds()
                      > timedelta(days=31).total_seconds()
                      for _start_date, _end_date in segments]):
            segments = [(from_datetime, to_datetime)
                        for from_datetime, to_datetime
                        in self.__get_datespan(start_date, end_date, interval)]
            interval += 1
        for _start_date, _end_date in segments:
            yield _start_date, _end_date

    def __get_datespan(self, startdate, enddate, interval):
        start_epoch = calendar.timegm(startdate.timetuple())
        end_epoch = calendar.timegm(enddate.timetuple())
        date_diff = end_epoch - start_epoch
        step = date_diff / interval
        delta = timedelta(seconds=step)
        currentdate = startdate
        while currentdate + delta <= enddate:
            todate = (currentdate + delta)
            yield currentdate, todate
            currentdate += delta

    @property
    def station(self):
        try:
            return self.__station
        except AttributeError:
            raise AttributeError('Must set station attribute.')

    @property
    def datetime(self):
        return self.__storage[self.station]['datetime']

    @property
    def values(self):
        return self.__storage[self.station]['values']

    @property
    def name(self):
        try:
            return self.__storage[self.station]['name']
        except KeyError:
            return ''

    @property
    def start_date(self):
        return self.__start_date

    @property
    def end_date(self):
        return self.__end_date

    @property
    def url(self):
        return "https://tidesandcurrents.noaa.gov/api/datagetter?"

    @property
    def datum(self):
        try:
            return self.__datum
        except AttributeError:
            return 'MSL'

    @property
    def units(self):
        try:
            return self.__units
        except AttributeError:
            return 'metric'

    @property
    def time_zone(self):
        try:
            return self.__time_zone
        except AttributeError:
            return 'gmt'

    @station.setter
    def station(self, station):
        assert station in self.__storage.keys()
        self.__station = station

    @start_date.setter
    def start_date(self, start_date):
        assert isinstance(start_date, datetime)
        self.__start_date = start_date

    @end_date.setter
    def end_date(self, end_date):
        assert isinstance(end_date, datetime)
        self.__end_date = end_date

    @datum.setter
    def datum(self, datum):
        assert datum \
            in ['MHHW', 'MHW', 'MTL', 'MSL', 'MLW', 'MLLW', 'NAVD88', 'STND']
        if datum == 'NAVD88':
            datum = 'NAVD'
        self.__datum = datum

    # def _call_REST(self):
    #     for station in self.stations:
    #         self._params['station'] = station
    #         response = requests.get(self._url, params=self._params)
    #         response.raise_for_status()
    #         data = json.loads(response.text)
    #         if "data" in data.keys():
    #             time = list()
    #             values=list()
    #             s=list()
    #             metadata=data['metadata']
    #             for datapoint in data['data']:
    #                 time.append(datetime.strptime(datapoint['t'], '%Y-%m-%d %H:%M'))
    #                 try:
    #                         val = float(datapoint['v'])
    #                 except:
    #                         val = np.nan
    #                 values.append(val)
    #                 try:
    #                         _s=float(datapoint['s'])
    #                 except:
    #                         _s=np.nan
    #                 s.append(_s)
    #             self[station] = { "time"     : np.asarray(time),
    #                                                 "zeta"     : np.ma.masked_invalid(values),
    #                                                 "s"        : np.ma.masked_invalid(s),
    #                                                 "metadata" : metadata,
    #                                                 "datum"    : self._params["datum"]}
