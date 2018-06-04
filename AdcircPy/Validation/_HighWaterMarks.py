import csv
import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from osgeo import ogr, osr
from AdcircPy import Validation

def from_csv(path):
    csvfile = open(path, 'r')
    f = csv.reader(csvfile)
    header = next(f)
    hwm_stations = dict()
    for line in f:
        station_id = line[header.index('hwm_id')]
        hwm_stations[station_id] = dict()
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
            del _HWM[station]
            cnt+=1
        if 'good' in type_list and self[station]['quality_id']==2:
            del _HWM[station]
            cnt+=1
        if 'fair' in type_list and self[station]['quality_id']==3:
            del _HWM[station]
            cnt+=1
        if 'poor' in type_list and self[station]['quality_id']==4:
            del _HWM[station]
            cnt+=1
        
        if 'riverine' in type_list and self[station]['environment_type'].lower()=='riverine':
            del _HWM[station]
            cnt+=1
    
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


def plot_HWM(self, HWM, extent=None, axes=None, vmin=None, vmax=None, title=None):
    if extent is None:
        extent = self.get_extent()
    idx, = np.where(np.logical_and(
                    np.logical_and(self.x>=extent[0], self.x<=extent[1]),
                    np.logical_and(self.y>=extent[2], self.y<=extent[3])))
    if vmin is None:
        vmin = np.min(self.values[idx])
    if vmax is None:
        vmax = np.max(self.values[idx])
    cmap = plt.get_cmap('jet')
    axes = self.make_plot(extent=extent, axes=axes, vmin=vmin, vmax=vmax, title=title)
    xyz = HWM.get_xyz()
    axes.scatter(xyz[:,0], xyz[:,1], c=xyz[:,2], vmin=vmin, vmax=vmax, cmap='jet')
    return axes