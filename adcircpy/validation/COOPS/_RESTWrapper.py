


class _RESTWrapper:

    def __init__(self, product, start_date, end_date, format='json',
                 units='metric', time_zone='gmt', datum='msl'):
        self._product = product

        self._init_params()
        self._call_REST()

    def __getitem__(self, key):
        return self._storage[key]

    def __iter__(self):
        return iter(self._storage)

    def __len__(self):
        return len(self._storage.keys())

    def _init_params(self):
        self._params={"format"       : "json",
                                     "units"       : "metric", 
                                     "time_zone"   : "gmt",
                                     "application" : "AdcircPy",
                                     "datum"       : "msl",
                                     "product"     : "water_level"}
        self._url = "https://tidesandcurrents.noaa.gov/api/datagetter?"
        self._params['begin_date'] = self.start_date.strftime('%Y%m%d %H:%M')
        self._params['end_date'] = self.end_date.strftime('%Y%m%d %H:%M')
        
    def _call_REST(self):
        for station in self.stations:
            self._params['station'] = station
            response = requests.get(self._url, params=self._params)
            response.raise_for_status()
            data = json.loads(response.text)
            if "data" in data.keys():
                time = list()
                values=list()
                s=list()
                metadata=data['metadata']
                for datapoint in data['data']:
                    time.append(datetime.strptime(datapoint['t'], '%Y-%m-%d %H:%M'))
                    try:
                            val = float(datapoint['v'])
                    except:
                            val = np.nan
                    values.append(val)
                    try:
                            _s=float(datapoint['s'])
                    except:
                            _s=np.nan
                    s.append(_s)
                self[station] = { "time"     : np.asarray(time),
                                                    "zeta"     : np.ma.masked_invalid(values),
                                                    "s"        : np.ma.masked_invalid(s),
                                                    "metadata" : metadata,
                                                    "datum"    : self._params["datum"]}
