from collections import OrderedDict
import csv
import json

import matplotlib.pyplot as plt
import requests


class HighWaterMarks(OrderedDict):
    url = 'https://stn.wim.usgs.gov/STNServices/HWMs/FilteredHWMs.json'
    params = {'EventType': 2,  # 2 for hurricane
              'EventStatus': 0}  # 0 for completed
    default_filter = {"riverine": True,
                      "non_still_water": True}

    def __init__(self, event_name, event_year,
                 filter_dict=default_filter, **station_data):
        super(HighWaterMarks, self).__init__(**station_data)
        self.event_name = event_name
        self.event_year = event_year
        self._filter(**filter_dict)

    @classmethod
    def from_event_name(cls, event_name, filter_dict=None):
        event_name, event_year, cls.params[
            'Event'] = cls._get_event_id_from_name(event_name)
        response = requests.get(cls.url, params=cls.params)
        response.raise_for_status()
        json_data = json.loads(response.text)
        hwm_stations = dict()
        for data in json_data:
            if 'elev_ft' in data.keys():
                hwm_stations[str(data['hwm_id'])] = data
        filter_dict = cls._init_filter_dict(filter_dict)
        return cls(event_name, event_year,
                   filter_dict=filter_dict, **hwm_stations)

    @classmethod
    def from_csv(cls, event_name, event_year, csvpath, filter_dict=None):
        filter_dict = cls._init_filter_dict(filter_dict)
        with open(csvpath, 'r') as f:
            csvfile = csv.reader(f)
            headers = next(csvfile)
            hwm_stations = dict()
            for i, line in enumerate(csvfile):
                station_id = str(i)
                hwm_stations[station_id] = dict()
                for j, column in enumerate(headers):
                    if column in ['longitude', 'latitude', 'elev_ft',
                                  'site_latitude', 'site_longitude']:
                        hwm_stations[station_id][column] = \
                            float(line[headers.index(column)])
                    elif column in ['hwm_quality_id']:
                        hwm_stations[station_id][column] = \
                            int(line[headers.index(column)])
                    else:
                        hwm_stations[station_id][column] = \
                            line[headers.index(column)]
        return cls(event_name, event_year,
                   filter_dict=filter_dict, **hwm_stations)

    def make_plot(self, axes=None, vmin=None, vmax=None,
                  extent=None, epsg=None, **kwargs):
        if axes is None:
            fig = plt.figure()
            axes = fig.add_subplot(111)
        for station in self.keys():
            axes.scatter(self[station]['lon'], self[station]['lat'],
                         c=self[station]['value'], vmin=vmin, vmax=vmax,
                         **kwargs)
        return axes

    @classmethod
    def _get_event_id_from_name(cls, event_name):
        """
        USGS is not standardizing event naming. Sometimes years are included,
        but sometimes they are ommitted. The order in which the names are in
        the response is also not standardized. Some workaround come into play
        in this algorithm in order to identify and categorize the dataset.
        USGS should standardize their REST server data.
        """
        response = requests.get(cls.url, params=cls.params)
        response.raise_for_status()
        json_data = json.loads(response.text)
        events = set()
        for item in json_data:
            event = item['eventName'].split()
            if len(event) == 1:
                eventName = event.pop()
                if eventName == 'Sandy':
                    eventYear = 2012
                elif eventName == 'Joaquin':
                    eventYear = 2015
                elif eventName == 'Hermine':
                    eventYear = 2016
                elif eventName == 'Irene':
                    eventYear = 2011
                elif eventName == 'Rita':
                    eventYear = 2005
                elif eventName == 'Wilma':
                    eventYear = 2005
                else:
                    eventYear = None
            elif len(event) == 2:
                eventName = event[0]
                eventYear = int(event[1])
            elif len(event) == 3:
                eventName = event.pop(0)
                for _ in event:
                    try:
                        eventYear = int(_)
                        break
                    except:
                        pass
            events.add((eventName, eventYear, int(item['event_id'])))
        events_dict = dict()
        for name, year, _id in events:
            events_dict[name.lower() + str(year)] = \
                {'name': name, 'year': year, 'id': _id}
        if event_name.lower() in events_dict.keys():
            return events_dict[event_name.lower()]['name'], \
                   events_dict[event_name.lower()]['year'], \
                   events_dict[event_name.lower()]['id']
        else:
            eventNames = [event.capitalize()
                          for event in list(events_dict.keys())]
            raise Exception('\nEvent name not Found!\n \
                            Valid event names are:\n{}'.format(eventNames))

    @classmethod
    def _init_filter_dict(cls, filter_dict):
        if filter_dict is None:
            filter_dict = cls.default_filter
        return filter_dict

    def _filter(self, excellent=False, good=False, fair=False, poor=False,
                riverine=False,
                non_still_water=False, keep_undefined=False):
        stations_to_delete = set()
        for station in self.keys():
            if 'hwm_quality_id' in self[station].keys():
                qid = self[station]['hwm_quality_id']
                if qid not in [1, 2, 3, 4]:
                    qid = None
            else:
                qid = None
            if qid is None:
                if keep_undefined == False:
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
                    stations_to_delete.add(station)
        for station in stations_to_delete:
            del self[station]
        self.filtered_count = len(stations_to_delete)
