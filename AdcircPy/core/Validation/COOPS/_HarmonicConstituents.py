from bs4 import BeautifulSoup
import json, codecs
import requests
from AdcircPy.core.Validation import COOPS
import os

_cachedir = os.getenv('LOCALAPPDATA')
if _cachedir is None:
    _cachedir = os.getenv('HOME')+'/.cache/AdcircPy'
else: 
    _cachedir += '/AdcircPy'
os.makedirs(_cachedir, exist_ok=True)

harmConstCache = _cachedir + '/harm_const.json'
if os.path.isfile(harmConstCache):
    harmConst = json.loads(open(harmConstCache, 'r').readlines()[0])
else:
    harmConst = dict()

rest_url = "https://tidesandcurrents.noaa.gov/api/datagetter?"
params = dict()
params['format'] = 'json'
params['unit'] = 'metric'
params['application'] = 'AdcircPy'
params['time_zone'] = 'gmt'


def _from_stations_list(stations):
  def __rebuild_cache():
    try:
      _rebuild_cache
      with open(harmConstCache, 'wb') as f:
        json.dump(harmConst, codecs.getwriter('utf-8')(f), ensure_ascii=False)
    except:
      pass
  # Note: Harmonic constituents not available through REST
  _stations = dict()
  for station in stations:
    if station not in harmConst.keys():
      _rebuild_cache = True
      url = 'https://tidesandcurrents.noaa.gov/harcon.html?id={}'.format(station)
      soup = BeautifulSoup(requests.get(url).text, 'html.parser')
      table = soup.find("table")
      # The first tr contains the field names.
      headings = [th.get_text().strip() for th in table.find("tr").find_all("th")]
      datasets = list()
      for row in table.find_all("tr")[1:]:
        datasets.append(dict(zip(headings, (td.get_text() for td in row.find_all("td")))))
      harmConst[station] = dict()
      for dataset in datasets:
        harmConst[station][dataset['Name']] = {
                    'amplitude'   : float(dataset['Amplitude'])/3.28084,
                    'phase'       : float(dataset['Phase']),
                    'speed'       : float(dataset['Speed']),
                    'description' : dataset['Description'],
                    'units'       : 'meters'}
    _stations[station]=harmConst[station]
  __rebuild_cache()
  return _stations

def _from_coops_id(station):
    if station not in harmConst.keys():
        pass
