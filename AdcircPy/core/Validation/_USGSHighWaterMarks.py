import argparse
from copy import deepcopy
from collections import defaultdict
import csv
import json
import matplotlib.pyplot as plt
from matplotlib.path import Path
import numpy as np
# from osgeo import ogr, osr
import psycopg2
import requests
from AdcircPy.Mesh import AdcircMesh
from AdcircPy.core import Validation

rest_url = 'https://stn.wim.usgs.gov/STNServices/HWMs/FilteredHWMs.json'
params = dict()
params['EventType'] = 2 # 2 for hurricane
params['EventStatus'] = 0 # for completed

def get_event_list():
    response = requests.get(rest_url, params=params)
    response.raise_for_status()
    json_data = json.loads(response.text)
    events=defaultdict(lambda : None)
    for item in json_data:
        events[item['eventName']]
    return [x.split()[0] for x in events.keys()]
     
def _get_event_id(eventName):
    response = requests.get(rest_url, params=params)
    response.raise_for_status()
    json_data = json.loads(response.text)
    events=defaultdict()
    for item in json_data:
        events[item['eventName'].split()[0].lower()]=item['event_id']
    events = dict(events)
    return events[eventName]

def from_event_name(eventName):
    params['Event'] = Validation.USGSHighWaterMarks._get_event_id(eventName.lower())
    response = requests.get(rest_url, params=params)
    response.raise_for_status()
    json_data = json.loads(response.text)
    hwm_stations = dict()
    for station_id, data in enumerate(json_data):
        station_id = str(station_id)
        hwm_stations[station_id] = dict()
        if 'elev_ft' in data.keys():
            for key in data.keys():
                hwm_stations[station_id][key] = data[key]
            hwm_stations[station_id]['elev_m'] = hwm_stations[station_id]['elev_ft'] / 3.28084
    return Validation.USGSHighWaterMarks(**hwm_stations)


def from_csv(path):
    csvfile = open(path, 'r')
    f = csv.reader(csvfile)
    header = next(f)
    hwm_stations = dict()
    for i, line in enumerate(f):
        station_id = str(i)
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
    csvfile.close()
    return Validation.HighWaterMarks(**hwm_stations)

def get_environments(self):
    environment = list()
    for station in self.keys():
        environment.append(self[station]['environment'])
    return environment
    
def filter(self, excellent=False, good=False, fair=False, poor=False, riverine=False, non_still_water=False, return_count=False, copy=True):
  if copy==True:
    _self = self
    self = deepcopy(_self)

  stations_to_delete = set()
  for station in self.keys():

    try:
        qid = self[station]['hwm_quality_id']
    except:
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
      if 'environment_type' in self[station].keys():
        if self[station]['environment_type'].lower() == 'riverine':
          stations_to_delete.add(station)

    if non_still_water == True and 'still_water' in self[station].keys():
      print('line 126 _USGSHighWaterMaks.py reports finding a still_water key.')
      print(self[station]['still_water'])
      stations_to_delete.add(station)        

  for station in stations_to_delete:
    del self[station]

  if return_count:
    return self, len(stations_to_delete)
  else:
    return self
        
def get_xy(self):
    lon = list()
    lat = list()
    for station in self.keys():
        lon.append(self[station]['site_longitude'])
        lat.append(self[station]['site_latitude'])    
    lon = np.asarray(lon).flatten()
    lat = np.asarray(lat).flatten()
    return np.vstack((lon, lat)).T

def get_values(self):
    values = list()
    for station in self.keys():
        values.append(self[station]['elev_m'])
    return np.array(values).flatten()

def get_xyz(self, units='meters'):
    lon = list()
    lat = list()
    values = list()
    for station in self.keys():
        lon.append(self[station]['site_longitude'])
        lat.append(self[station]['site_latitude'])  
        if units=='meters':
            values.append(self[station]['elev_m']) 
        elif units=='feet':
            values.append(self[station]['elev_ft']) 
    return np.vstack((lon,lat,values)).T

def get_counties(self):
    counties=list()
    for key in self.keys():
        counties.append(self[key]['county'])
    return set(counties)

def export_shapefile(self, path, layer_name='High Water Marks', epsg=4326):
    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.CreateDataSource(path)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    layer = data_source.CreateLayer(layer_name, srs, ogr.wkbPoint)
    layer.CreateField(ogr.FieldDefn("Latitude", ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn("Longitude", ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn("Elevation", ogr.OFTReal))
    for station in self.keys():
        feature = ogr.Feature(layer.GetLayerDefn())
        feature.SetField("Latitude", self[station]['lat'])
        feature.SetField("Longitude", self[station]['lon'])
        feature.SetField("Elevation", self[station]['value'])
        wkt = "POINT(%f %f)" %  (self[station]['lon'] , self[station]['lat'])
        point = ogr.CreateGeometryFromWkt(wkt)
        feature.SetGeometry(point)
        layer.CreateFeature(feature)
        feature = None

def export_to_PostGIS(self, dbname, **kwargs):
    user        = kwargs.pop('user', 'postgres')
    password    = kwargs.pop('password', True)
    host        = kwargs.pop('host', 'localhost')
    port        = kwargs.pop('port', 5432)
    overwrite   = kwargs.pop('overwrite', False)
    schema      = kwargs.pop('schema', 'ADCIRC')
    table       = kwargs.pop('table', 'high_water_marks')
    if password==True:
        password = getpass.getpass('Password: ')
    con = psycopg2.connect(dbname=dbname, user=user, password=password, host=host, port=port)
    cur  = con.cursor()
    cur.execute('CREATE SCHEMA IF NOT EXISTS {};'.format(schema))
    con.commit()
    cur.execute("SELECT to_regclass('{}.{}');".format(schema, table))
    table_existance = cur.fetchall()[0][0]
    if table_existance is None or overwrite==True:
        if table_existance is not None and overwrite==True:
            cur.execute('DROP TABLE {}.{};'.format(schema, table))
            con.commit()
        cur.execute('CREATE TABLE IF NOT EXISTS {}.{} (id SERIAL PRIMARY KEY, geom geometry(Point, {}), value REAL);'.format(schema, table, self.epsg))
        for key in self.keys():
            geom = "ST_GeomFromText('Point({:f} {:f})', {})".format(self[key]['lon'], self[key]['lat'], self.epsg)
            cur.execute("INSERT INTO {}.{} (geom, value) VALUES ({}, {});".format(schema, table, geom, self[key]['value']))
        con.commit()
        cur.execute("CREATE INDEX sidx_{}_geom ON {}.{} USING GIST (geom);".format(table, schema, table))
        con.commit()
    elif table_existance is not None and overwrite==False:
        raise Exception("Table {}.{} exists in database. Manually delete or use overwrite=True".format(schema, table))


def clip_from_shapefile(self, path, return_count=False):
    shapefile = ogr.Open(path)
    layer = shapefile.GetLayer()
    polygons=list()
    for area in layer: 
        area_shape = area.GetGeometryRef() 
        area_polygon = area_shape.GetGeometryRef(0) 
        no_of_polygon_vertices = area_polygon.GetPointCount() 
        vertices = list()
        codes = [Path.MOVETO]
        for vertex in range(no_of_polygon_vertices):
            lon, lat, z = area_polygon.GetPoint(vertex)
            vertices.append((lon, lat))
            codes.append(Path.LINETO)
        vertices.append(vertices[0])
        codes[-1] = Path.CLOSEPOLY
        polygons.append(Path(vertices, codes))
    to_remove=list()
    for station in self.keys():
        coords = (self[station]['lon'], self[station]['lat'])
        for polygon in polygons:
            if polygon.contains_point(coords):
                to_remove.append(station)
    _hwm = copy.deepcopy(self)
    for station in self.keys():
        if station in to_remove:
            del _hwm[station]
    if return_count==True:
        return _hwm, len(to_remove)
    else:
        return _hwm


def make_plot(self, axes=None, vmin=None, vmax=None, extent=None, epsg=None, **kwargs):
    if axes is None:                
        fig = plt.figure()
        axes  = fig.add_subplot(111)
    for station in self.keys():
        axes.scatter(self[station]['lon'], self[station]['lat'], c=self[station]['value'], vmin=vmin, vmax=vmax, **kwargs)
    return axes

def get_from_extent(self, extent, epsg):
    hwm = copy.deepcopy(self)
    path = Path([(extent[0], extent[2]),
                 (extent[1], extent[2]),
                 (extent[1], extent[3]),
                 (extent[0], extent[3]),
                 (extent[0], extent[2])], closed=True)

    if epsg != 4326:
        target_proj = pyproj.Proj(init="epsg:{}".format(epsg))
    
    for station in list(hwm.keys()):
        x = hwm[station]['site_longitude']
        y = hwm[station]['site_latitude']
        if epsg != 4326:
            x, y = target_proj(x, y, inverse=True)
        if path.contains_point((x, y))==False:
            hwm.pop(station)
    return hwm

