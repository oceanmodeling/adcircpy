import argparse
import csv
import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path

import psycopg2
from osgeo import ogr, osr
from AdcircPy import Validation

def from_csv(path):
    csvfile = open(path, 'r')
    f = csv.reader(csvfile)
    header = next(f)
    hwm_stations = dict()
    for i, line in enumerate(f):
        station_id = str(i)
        hwm_stations[station_id] = dict()
        hwm_stations[station_id]['hwm_id'] = line[header.index('hwm_id')]
        hwm_stations[station_id]['lat']=float(line[header.index('site_latitude')])
        hwm_stations[station_id]['lon']=float(line[header.index('site_longitude')])
        hwm_stations[station_id]['state'] = line[header.index('stateName')]
        hwm_stations[station_id]['county'] = line[header.index('countyName')]
        hwm_stations[station_id]['location_name'] = line[header.index('countyName')] +', '+line[header.index('stateName')]
        hwm_stations[station_id]['site_description'] = line[header.index('hwm_locationdescription')]
        hwm_stations[station_id]['quality'] = line[header.index('hwmQualityName')]
        hwm_stations[station_id]['quality_id'] = int(line[header.index('hwm_quality_id')])
        hwm_stations[station_id]['value'] = float(line[header.index('elev_ft')])
        hwm_stations[station_id]['units'] = 'feet'
        hwm_stations[station_id]['vertical_datum'] = line[header.index('verticalDatumName')]
        hwm_stations[station_id]['horizontal_datum'] = line[header.index('horizontalDatumName')]
        hwm_stations[station_id]['environment_type'] = line[header.index('hwm_environment')].lower()
        try:
            still_water = int(line[header.index('stillwater')])
            if still_water == 1:
                still_water = True
            else:
                still_water = False
        except:
            still_water = None
        hwm_stations[station_id]['still_water'] = still_water
    csvfile.close()
    return Validation.HighWaterMarks(**hwm_stations)

def _Validate(self, Mesh):
    raise NotImplementedError

def from_event_id(id):
	raise NotImplementedError("Coming soon!")


def get_environments(self):
    environment = list()
    for station in self.keys():
        environment.append(self[station]['environment'])
    return environment
    
def remove(self, type_list, return_count=False):
    _HWM = copy.deepcopy(self)
    try:
        type_list = type_list.split()
    except:
        pass
    type_list = [x.lower() for x in type_list]
    cnt=0
    for station in self.keys():
        if 'excellent' in type_list and self[station]['quality_id']==1:
            try:
                del _HWM[station]
                cnt+=1
            except:
                pass
        if 'good' in type_list and self[station]['quality_id']==2:
            try:
                del _HWM[station]
                cnt+=1
            except:
                pass
        if 'fair' in type_list and self[station]['quality_id']==3:
            try:
                del _HWM[station]
                cnt+=1
            except:
                pass
        if 'poor' in type_list and self[station]['quality_id']==4:
            try:   
                del _HWM[station]
                cnt+=1
            except:
                pass
        if 'riverine' in type_list and self[station]['environment_type'].lower()=='riverine':
            try:
                del _HWM[station]
                cnt+=1
            except:
                pass
    if return_count:
        return _HWM, cnt
    else:
        return _HWM
        
def still_water_only(self, return_count=False, keep_undefined=False):
    _HWM = copy.deepcopy(self)
    cnt=0
    for station in self.keys():
        if self[station]['still_water']==False:
            del _HWM[station]
            cnt+=1
        if keep_undefined==False:
            if self[station]['still_water']==None:
                del _HWM[station]
                cnt+=1
    if return_count:
        return _HWM, cnt
    else:
        return _HWM

def convert_to_meters(self):
    for station in self.keys():
        if self[station]['units'] == 'feet':
            self[station]['value'] /= 3.28084
            self[station]['units'] == 'meters'

def get_coordinates(self):
    lon = list()
    lat = list()
    for station in self.keys():
        lon.append(self[station]['lon'])
        lat.append(self[station]['lat'])    
    lon = np.asarray(lon).flatten()
    lat = np.asarray(lat).flatten()
    return lon, lat

def get_values(self):
    values = list()
    for station in self.keys():
        values.append(self[station]['value'])
    return np.array(values).flatten()

def get_xyz(self):
    lon, lat = self.get_coordinates()
    values = self.get_values()
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
        x = hwm[station]['lon']
        y = hwm[station]['lat']
        if epsg != 4326:
            x, y = target_proj(x, y, inverse=True)
        if path.contains_point((x, y))==False:
            hwm.pop(station)
    return hwm

def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('wavemodel', choices=['surfbeat', 'nonh'], help='Wave model type: surfbeat or non-hydrostatic')
    parser.add_argument('grid_id', type=str, help='Grid grid id')
    parser.add_argument('CAT', choices=[1,2,3,4,5], type=int, help='Hurricane Category')
    parser.add_argument('--SLR', '-s', choices=[0., 0.5, 1.], dest='SLR', type=float, default=0., help='Sea Level Rise Value')
    parser.add_argument('--wavedir', '-wd', type=float, help='Optional wave direction in Nautical Coordinates. Wavemaker aligned to offshore boundary assumed when not given.')
    parser.add_argument('--tma', choices=[0, 1], default=1, type=int, help='TMA vs Jonswap spectrum. 0 for Jonswap, 1 for TMA (default)')
    parser.add_argument('--dx', '-dx', default=None, type=float, help='Zonal resolution in meters. Defaults to 3 meters for nonh and 10 meters for surfbeat.')
    parser.add_argument('--dy', '-dy', default=None, type=float, help='Meridional resolution in meters. Defaults to 3 meters for nonh and 10 meters for surfbeat.')
    parser.add_argument('--epsg', '-e', default=32620, type=int, help='EPSG code')
    parser.add_argument('--spread', '-sp', default=15., type=float, help='Wave angular spreading in degrees.')
    parser.add_argument('--ncpu', '-n', default=600, type=int, help='Total number of CPUs to use.')
    parser.add_argument('--filter', '-f', default=2, type=int, dest='sigma', help='Apply gaussian filter to bathy.')
    parser.add_argument('--thetanaut', choices=[0,1], type=int, default=1, help='Advanced option. Wave direction given in nautical convention (1) or angle wrt to grid x-axis (0).')
    parser.add_argument('--thetamin', type=float,  help='Advanced option. Thetamin directional bin for surfbeat model')
    parser.add_argument('--thetamax', type=float, help='Advanced option. Thetamax directional bin for surfbeat model')
    parser.add_argument('--dtheta', default=10., type=float, help='Advanced option. Thetamax for directional bin for surfbeat model')
    parser.add_argument('--ts-length', '-rt', type=int, dest='rt', help='Advanced option. Time series record length rt.')
    self.args = parser.parse_args()


def _main():
    self._parse_args()

