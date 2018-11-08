from collections import OrderedDict, defaultdict
import requests
import json

class HighWaterMarks(OrderedDict):
  url = 'https://stn.wim.usgs.gov/STNServices/HWMs/FilteredHWMs.json'
  params =  {'EventType'   : 2, # 2 for hurricane
             'EventStatus' : 0} # 0 for completed
  
  def __call__(self, **station_data):
    super(HighWaterMarks, self).__init__(**station_data)

  @classmethod
  def from_event_name(cls, event_name, target_datum, filter_dict=None):
    cls.params['Event'] = cls.get_event_id_from_name(event_name)
    if filter_dict is None:
      filter_dict  = { "riverine" : True,
                       "non_still_water" : True }
    response = requests.get(cls.url, params=cls.params)
    response.raise_for_status()
    json_data= json.loads(response.text)
    hwm_stations=dict()
    for data in json_data:
      if 'elev_ft' in data.keys():
        hwm_stations[str(data['hwm_id'])]=data
    cls = cls(**hwm_stations)
    cls.filter(**filter_dict)
    return cls

  @classmethod
  def from_csv(cls, path):
    csvfile = open(path, 'r')
    with csv.reader(csvfile) as f: 
      header = next(f)
      hwm_stations = dict()
      for i, line in enumerate(f):
        # station_id = str(i)
        hwm_stations[station_id] = dict()
        hwm_stations[station_id]['site_latitude']=float(line[header.index('site_latitude')])
        hwm_stations[station_id]['site_longitude']=float(line[header.index('site_longitude')])
        hwm_stations[station_id]['stateName'] = line[header.index('stateName')]
        hwm_stations[station_id]['countyName'] = line[header.index('countyName')]
        hwm_stations[station_id]['hwm_locationdescription'] = line[header.index('hwm_locationdescription')]
        hwm_stations[station_id]['hwmQualityName'] = line[header.index('hwmQualityName')]
        hwm_stations[station_id]['hwm_quality_id'] = int(line[header.index('hwm_quality_id')])
        hwm_stations[station_id]['elev_ft'] = float(line[header.index('elev_ft')])
        hwm_stations[station_id]['verticalDatumName'] = line[header.index('verticalDatumName')]
        hwm_stations[station_id]['horizontalDatumName'] = line[header.index('horizontalDatumName')]
        hwm_stations[station_id]['hwm_environment'] = line[header.index('hwm_environment')].lower()
        hwm_stations[station_id]['elev_m'] = hwm_stations[station_id]['elev_ft'] / 3.28084
    return cls(**hwm_stations)


  @classmethod
  def get_event_id_from_name(cls, eventName):
    response = requests.get(cls.url, params=cls.params)
    response.raise_for_status()
    json_data = json.loads(response.text)
    events=defaultdict()
    for item in json_data:
        events[item['eventName'].split()[0].lower()]=item['event_id']
    events = dict(events)
    if eventName.lower() in events.keys():
        return events[eventName.lower()]
    else:
        raise Exception('eventName not Found! Valid event names are: {}'.format(list(events.keys())))

  def filter(self, excellent=False, good=False, fair=False, poor=False, riverine=False, non_still_water=False, keep_undefined=False):
    stations_to_delete = set()
    for station in self.keys():
      if 'hwm_quality_id' in self[station].keys():
        qid = self[station]['hwm_quality_id']
        if qid not in [1,2,3,4]:
          qid=None
      else:
        qid=None
      if qid is None:         
        if keep_undefined==False:
          stations_to_delete.add(station)
      if excellent == True and qid == 1:
        stations_to_delete.add(station)
      if good == True and qid == 2:
        stations_to_delete.add(station)
      if fair == True and qid == 3:
        stations_to_delete.add(station)
      if poor == True and qid == 4:
        stations_to_delete.add(station)
      if riverine == True:
        if 'hwm_environment' in self[station].keys():
          if self[station]['hwm_environment'].lower() == 'riverine':
            stations_to_delete.add(station)
      if non_still_water == True:
          if 'still_water' in self[station].keys():
              print('_USGSHighWaterMaks.py reports finding a still_water key on filter() (report this to the devs! jreniel@gmail.com)')
              print(self[station]['still_water'])
              stations_to_delete.add(station)        
    for station in stations_to_delete:
      del self[station]
    self.filtered_count=len(stations_to_delete)
 
  def make_plot(self, axes=None, vmin=None, vmax=None, extent=None, epsg=None,**kwargs):
    if axes is None:                
        fig = plt.figure()
        axes  = fig.add_subplot(111)
    for station in self.keys():
        axes.scatter(self[station]['lon'], self[station]['lat'], c=self[station]['value'], vmin=vmin, vmax=vmax, **kwargs)
    return axes