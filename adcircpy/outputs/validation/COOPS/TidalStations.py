from collections.abc import Mapping
import numpy as np
from datetime import datetime
import json
import requests
from AdcircPy.Outputs._StationTimeseries import _StationTimeseries


class TidalStations(Mapping):

    def __init__(self, _format='json', units='metric', time_zone='gmt',
                 datum='msl'):
        self._storage = dict()
        self._url = "https://tidesandcurrents.noaa.gov/api/datagetter?"
        self._params = {}
        self._product = 'water_level'
        self._format = _format
        self._units = units
        self._time_zone = time_zone
        self._datum = datum

    def __getitem__(self, key):
        return self._storage[key]

    def __iter__(self):
        return iter(self._storage)

    def __len__(self):
        return len(self._storage.keys())

    def fetch(self, station_id, start_date, end_date):
        station_id = str(station_id)
        assert isinstance(start_date, datetime)
        assert isinstance(end_date, datetime)
        params = self.params
        params['station'] = station_id
        # dt = end_date - start_date
        # if dt.total_seconds() > 30.*24.*60.*60.:
        #     _end_date = _start_date + deltatime()
        #     while dt.total_seconds() > 30.*24.*60.*60.:

        # else:
        params['begin_date'] = start_date.strftime('%Y%m%d')
        params['end_date'] = end_date.strftime('%Y%m%d')
        response = requests.get(self._url, params=self._params)
        response.raise_for_status()
        json_data = json.loads(response.text)
        if 'error' in json_data.keys():
            pass
        else:
            x = json_data['metadata']['lon']
            y = json_data['metadata']['lat']
            name = json_data['metadata']['name']
            time = list()
            values = list()
            for data in json_data['data']:
                date = datetime.strptime(data['t'], '%Y-%m-%d %H:%M')
                if date >= start_date and date <= end_date:
                    time.append(date)
                    try:
                        values.append(float(data['v']))
                    except ValueError:
                        values.append(np.nan)
            self.add_station(station_id, x, y, values, time, name)

    def add_station(self, station_id, x, y, values, time, name):
        self._storage[station_id] = _StationTimeseries(
                                                    x, y, values, time, name)

    @property
    def url(self):
        return self._url

    @property
    def params(self):
        return self._params

    @property
    def _storage(self):
        return self.__storage

    @property
    def _url(self):
        return self.__url

    @property
    def _params(self):
        return self.__params

    @property
    def _product(self):
        return self.__product

    @property
    def _format(self):
        return self.__format

    @property
    def _units(self):
        return self.__units

    @property
    def _time_zone(self):
        return self.__time_zone

    @property
    def _datum(self):
        return self.__datum

    @_storage.setter
    def _storage(self, storage):
        self.__storage = dict()

    @_url.setter
    def _url(self, url):
        self.__url = url

    @_params.setter
    def _params(self, params):
        self.__params = {}

    @_product.setter
    def _product(self, product):
        self.params['product'] = product

    @_format.setter
    def _format(self, _format):
        self.params['format'] = _format

    @_units.setter
    def _units(self, units):
        self.params['units'] = units

    @_time_zone.setter
    def _time_zone(self, time_zone):
        self.params['time_zone'] = time_zone

    @_datum.setter
    def _datum(self, datum):
        self.params['datum'] = datum


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

