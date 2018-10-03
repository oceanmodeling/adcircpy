from copy import deepcopy
from collections import defaultdict
import csv
import tqdm
import json
import matplotlib.pyplot as plt
from matplotlib.path import Path
import numpy as np
import requests
import warnings
from AdcircPy.Datum import VDatum


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
  cls._filter(**filter_dict)
  cls._set_vertical_datum(target_datum)
  return cls


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
  cls(**hwm_stations)


def get_environments(self):
    environment = list()
    for station in self.keys():
        environment.append(self[station]['environment'])
    return environment
    
def _filter(self, excellent=False, good=False, fair=False, poor=False, riverine=False, non_still_water=False, keep_undefined=False):
  # print(self.keys())
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
   

def _set_vertical_datum(self, target_datum, _tqdm=False):
  data_by_hdatum=defaultdict(list)
  for station in self.keys():
    s_h_datum = self[station]['horizontalDatumName']
    if 'WGS84' in s_h_datum:
      # See rule: WGS84_TRANSIT   WGS84(transit) - use NAD83 (see NGS's HTDP)
      # 'WGS84' --> 'NAD83' for practical purposes.
      s_h_datum='NAD83'
    data_by_hdatum[s_h_datum].append(self[station])

  for s_h_datum in data_by_hdatum.keys():
    data_by_vdatum=defaultdict(list)
    for dataset in data_by_hdatum[s_h_datum]:
      data_by_vdatum[dataset['verticalDatumName']].append(dataset)
    
    for s_v_datum in data_by_vdatum.keys():
      if len(self.keys())>200:
        warnings.warn("Because data is being streamed from a REST service, this process may hang unexpectedly for large datasets. \
                       Use _tqdm=True keyword argument to see the progress bar.")
      if _tqdm==True:
        print('Converting {}:{} data to NAD83:LMSL, please wait...'.format(s_h_datum, s_v_datum))
        iterable=tqdm.tqdm()
      else:
        iterable=data_by_vdatum[s_v_datum]

      for data in iterable:
        _x = self[str(data['hwm_id'])]['longitude']
        _y = self[str(data['hwm_id'])]['latitude']
        _z = self[str(data['hwm_id'])]['elev_ft']
        try:
          if s_v_datum.strip('\n ') != target_datum.strip('\n '):
          # Some combinations like NAD83:NGVD29 to NAD83:LMSL will return
          # "Source Horizontal Frame should be NAD27" or similar message.
            value = VDatum(_x, _y, _z, s_h_datum, s_v_datum, 'us_ft', target_datum)
          else:
            value=_z/3.28084
          if value != -999999.:
            self[str(data['hwm_id'])]['elev_m'] = value
            self[str(data['hwm_id'])]['verticalDatumName'] = target_datum
          elif value==-999999. or value==-99999./3.28084:
            del self[str(data['hwm_id'])]
            self.filtered_count+=1
        except Exception:
          tqdm.tqdm.write("Datapoint with hwm_id={} could not be processed.".format(str(data['hwm_id'])))
          del self[str(data['hwm_id'])]
          self.filtered_count+=1


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
    counties=set()
    for key in self.keys():
        counties.add(self[key]['county'])
    return counties

def make_plot(self, axes=None, vmin=None, vmax=None, extent=None, epsg=None, **kwargs):
    if axes is None:                
        fig = plt.figure()
        axes  = fig.add_subplot(111)
    for station in self.keys():
        axes.scatter(self[station]['lon'], self[station]['lat'], c=self[station]['value'], vmin=vmin, vmax=vmax, **kwargs)
    return axes

def _clip_to_extent(self, extent, epsg):
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

def get_event_id_from_name(cls, eventName):
  # cls.__init__(cls)
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

# def export_shapefile(self, path, layer_name='High Water Marks', epsg=4326):
#     driver = ogr.GetDriverByName("ESRI Shapefile")
#     data_source = driver.CreateDataSource(path)
#     srs = osr.SpatialReference()
#     srs.ImportFromEPSG(epsg)
#     layer = data_source.CreateLayer(layer_name, srs, ogr.wkbPoint)
#     layer.CreateField(ogr.FieldDefn("Latitude", ogr.OFTReal))
#     layer.CreateField(ogr.FieldDefn("Longitude", ogr.OFTReal))
#     layer.CreateField(ogr.FieldDefn("Elevation", ogr.OFTReal))
#     for station in self.keys():
#         feature = ogr.Feature(layer.GetLayerDefn())
#         feature.SetField("Latitude", self[station]['lat'])
#         feature.SetField("Longitude", self[station]['lon'])
#         feature.SetField("Elevation", self[station]['value'])
#         wkt = "POINT(%f %f)" %  (self[station]['lon'] , self[station]['lat'])
#         point = ogr.CreateGeometryFromWkt(wkt)
#         feature.SetGeometry(point)
#         layer.CreateFeature(feature)
#         feature = None

# def export_to_PostGIS(self, dbname, **kwargs):
#     user        = kwargs.pop('user', 'postgres')
#     password    = kwargs.pop('password', True)
#     host        = kwargs.pop('host', 'localhost')
#     port        = kwargs.pop('port', 5432)
#     overwrite   = kwargs.pop('overwrite', False)
#     schema      = kwargs.pop('schema', 'ADCIRC')
#     table       = kwargs.pop('table', 'high_water_marks')
#     if password==True:
#         password = getpass.getpass('Password: ')
#     con = psycopg2.connect(dbname=dbname, user=user, password=password, host=host, port=port)
#     cur  = con.cursor()
#     cur.execute('CREATE SCHEMA IF NOT EXISTS {};'.format(schema))
#     con.commit()
#     cur.execute("SELECT to_regclass('{}.{}');".format(schema, table))
#     table_existance = cur.fetchall()[0][0]
#     if table_existance is None or overwrite==True:
#         if table_existance is not None and overwrite==True:
#             cur.execute('DROP TABLE {}.{};'.format(schema, table))
#             con.commit()
#         cur.execute('CREATE TABLE IF NOT EXISTS {}.{} (id SERIAL PRIMARY KEY, geom geometry(Point, {}), value REAL);'.format(schema, table, self.epsg))
#         for key in self.keys():
#             geom = "ST_GeomFromText('Point({:f} {:f})', {})".format(self[key]['lon'], self[key]['lat'], self.epsg)
#             cur.execute("INSERT INTO {}.{} (geom, value) VALUES ({}, {});".format(schema, table, geom, self[key]['value']))
#         con.commit()
#         cur.execute("CREATE INDEX sidx_{}_geom ON {}.{} USING GIST (geom);".format(table, schema, table))
#         con.commit()
#     elif table_existance is not None and overwrite==False:
#         raise Exception("Table {}.{} exists in database. Manually delete or use overwrite=True".format(schema, table))


# def clip_from_shapefile(self, path, return_count=False):
#     shapefile = ogr.Open(path)
#     layer = shapefile.GetLayer()
#     polygons=list()
#     for area in layer: 
#         area_shape = area.GetGeometryRef() 
#         area_polygon = area_shape.GetGeometryRef(0) 
#         no_of_polygon_vertices = area_polygon.GetPointCount() 
#         vertices = list()
#         codes = [Path.MOVETO]
#         for vertex in range(no_of_polygon_vertices):
#             lon, lat, z = area_polygon.GetPoint(vertex)
#             vertices.append((lon, lat))
#             codes.append(Path.LINETO)
#         vertices.append(vertices[0])
#         codes[-1] = Path.CLOSEPOLY
#         polygons.append(Path(vertices, codes))
#     to_remove=list()
#     for station in self.keys():
#         coords = (self[station]['lon'], self[station]['lat'])
#         for polygon in polygons:
#             if polygon.contains_point(coords):
#                 to_remove.append(station)
#     _hwm = copy.deepcopy(self)
#     for station in self.keys():
#         if station in to_remove:
#             del _hwm[station]
#     if return_count==True:
#         return _hwm, len(to_remove)
#     else:
#         return _hwm