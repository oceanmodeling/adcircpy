from datetime import datetime, timedelta
import json
import numpy as np
import requests
from AdcircPy import COOPS

class _REST(object):
    _url = "https://tidesandcurrents.noaa.gov/api/datagetter?"
    _params={"format"       : "json",
             "units"         : "metric", 
             "time_zone"    : "gmt",
             "application" : "AdcircPy",
             "datum"       : "msl"}

    def __init__(self, **kwargs):
        self._products = ['air_gap',
                          'air_pressure',
                          'air_temperature',
                          'conductivity',
                          'currents',
                          'currents_survey',
                          'daily_mean',
                          'datums',
                          'high_low',
                          'hourly_height',
                          'humidity',
                          'monthly_mean',
                          'one_minute_water_level',
                          'predictions',
                          'salinity',
                          'visibility',
                          'water_level',
                          'water_temperature',
                          'wind']
        self._format = 'json'
        self._unit = 'metric'
        self._time_zone = 'gmt'
        self._application = 'AdcircPy'
        # self.params = {
        #             }

def _from_station_list(station_list, start_date, end_date):
    _REST._params['product'] = 'water_level'
    _REST._params['begin_date'] = start_date.strftime('%Y%m%d %H:%M')
    _REST._params['end_date'] = end_date.strftime('%Y%m%d %H:%M')
    _stations = dict()
    for station in station_list:
        _REST._params['station'] = station
        response = requests.get(_REST._url, params=_REST._params)
        response.raise_for_status()
        data = json.loads(response.text)

        if "data" in data.keys():
            time = list()
            values=list()
            s=list()
            metadata=data['metadata']
            for _datapoint in data['data']:
                time.append(datetime.strptime(_datapoint['t'], '%Y-%m-%d %H:%M'))
                try:
                    _val = float(_datapoint['v'])
                except:
                    _val = np.nan
                values.append(_val)
                try:
                    _s=float(_datapoint['s'])
                except:
                    _s=np.nan
                s.append(_s)
        else:
            time=[]
            values=[]
            s=[]
            metadata={'name':'', 'lon':'', 'lat':''}
        _stations[station] = { "time"      : np.asarray(time),
                                "zeta"     : np.ma.masked_invalid(values),
                                "s"        : np.ma.masked_invalid(s),
                                "metadata" : metadata,
                                "datum"    : _REST._params["datum"]}
    return COOPS.TidalStations(**_stations)