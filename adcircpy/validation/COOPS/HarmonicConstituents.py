from bs4 import BeautifulSoup
import json
import codecs
import requests
import os
from AdcircPy import utils


class HarmonicConstituents(dict):

    _rebuild = False

    def __init__(self, stations):
        super(HarmonicConstituents, self).__init__()
        self._stations = stations
        self._url = 'https://tidesandcurrents.noaa.gov/harcon.html?id='
        self._init_cache()
        self._init_stations()

    def _init_cache(self):
        self._cachedir = utils.get_cache_dir()
        self._harmConstCache = self._cachedir + '/harm_const.json'
        if os.path.isfile(self._harmConstCache):
            self._cache = json.loads(open(self._harmConstCache, 'r')
                                     .readlines()[0])
        else:
            self._cache = dict()

    def _init_stations(self):
        for station in self._stations:
            if station not in list(self._cache.keys()):
                self._add_station_to_cache(station)
            self[station] = self._cache[station]

    def _add_station_to_cache(self, station):
        # Parse from html, not available through rest.
        url = self._url+station
        soup = BeautifulSoup(requests.get(url).text, 'html.parser')
        table = soup.find("table")
        if table is not None:
            headings = [th.get_text().strip() for th in table.find("tr")
                        .find_all("th")]
            datasets = list()
            for row in table.find_all("tr")[1:]:
                datasets.append(dict(zip(headings, (td.get_text() for td in
                                                    row.find_all("td")))))
            for dataset in datasets:
                if station not in self._cache.keys():
                    self._cache[station] = dict()
                if float(dataset['Amplitude']) != 0.:
                    self._cache[station][dataset['Name']] = {
                            'amplitude': float(dataset['Amplitude'])/3.28084,
                            'phase': float(dataset['Phase']),
                            'speed': float(dataset['Speed']),
                            'description': dataset['Description'],
                            'units': 'meters'}
        else:
            self._cache[station] is None
        self._rebuild = True

    def __del__(self):
        if self._rebuild is True:
            with open(self._harmConstCache, 'wb') as f:
                json.dump(self._cache, codecs.getwriter('utf-8')(f),
                          ensure_ascii=False)
