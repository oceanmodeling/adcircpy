from __future__ import absolute_import, division, print_function
import os
import copy
import numpy as np
import matplotlib
from scipy.interpolate import RectBivariateSpline, griddata
from scipy.ndimage.filters import generic_filter, gaussian_filter
from haversine import haversine
from osgeo import gdal, osr
import pyproj
import utm
import psycopg2
import getpass
from AdcircPy import demtools
gdal.UseExceptions()

def read_tile (path, force_epsg=None):
    """
    Reads any DEM readable by GDAL and sets the attributes described in the 
    following list:

        self.x
        self.y
        self.values 
        self.epsg
        self.geoTransform

        geoTransform definition:
            self.geoTransform[0] /* top left x */
            self.geoTransform[1] /* w-e pixel resolution */
            self.geoTransform[2] /* 0 */
            self.geoTransform[3] /* top left y */
            self.geoTransform[4] /* 0 */
            self.geoTransform[5] /* n-s pixel resolution (negative value) */
    """
    geo  = gdal.Open(path)
    wkt = geo.GetProjection()
    inproj = osr.SpatialReference(wkt=wkt)
    vertical_datum = inproj.GetAttrValue('datum')
    if force_epsg!=None:
        self.epsg=epsg
    else:
        if inproj.IsGeographic() == 1:
            epsg = 4326
        else:
            epsg = inproj.GetAttrValue('PROJCS|AUTHORITY', 1)
            # Hack to force identification of WGS84 on certain files.
            if epsg is None and geo.GetGeoTransform()[2]<10**-4:
                epsg=4326
            elif epsg is None:
                raise Exception("Could not auto identify tile EPSG. You can use the force_epg kwarg.\nwkt is {}".format(inproj))
    
    geoTransform  = geo.GetGeoTransform()
    x    = np.linspace(geoTransform[0], 
                        geoTransform[0] + geo.RasterXSize * geoTransform[1],
                        num = geo.RasterXSize)
    y    = np.linspace(geoTransform[3] + geo.RasterYSize * geoTransform[5],
                        geoTransform[3],
                        num = geo.RasterYSize)
    values = geo.ReadAsArray()
    values = np.ma.masked_equal(values, geo.GetRasterBand(1).GetNoDataValue())
    return demtools.DEM(x, y, values, geoTransform, epsg, vertical_datum)


def circular_filter(self, radius):
    """
    Input radius in meters.
    """
    yres=self.get_dy(meters=True)
    xres=self.get_dx(meters=True)
    num = 0.
    while num < radius: num = num + xres
    x_radius = num - xres
    num = 0.
    while num < radius: num = num + yres
    y_radius = num - yres
    x = np.arange(-x_radius, x_radius+xres, xres)
    y = np.arange(-y_radius, y_radius+yres, yres)
    y,x = np.ogrid[-radius:radius+yres:yres, -radius:radius+xres:xres]
    mask = x**2 + y**2 <= radius**2
    kernel = np.zeros((mask.shape[0], mask.shape[1]))
    kernel[mask] = 1
    self.values = generic_filter(self.values, np.mean, footprint=kernel)

def circular_pixel_filter(self, radius):
    """
    Input radius in meters.
    """
    ncells = int(radius // resolution_in_meters(self))
    kernel = np.zeros((2*ncells+1, 2*ncells+1))
    y,x = np.ogrid[-ncells:ncells+1, -ncells:ncells+1]
    mask = x**2 + y**2 <= ncells**2
    kernel[mask] = 1
    self.values = generic_filter(self.values, np.mean, footprint=kernel)

def export_to_file(self, path, driver='Gtiff'):  
    gdal.UseExceptions()
    driver = gdal.GetDriverByName(driver)
    ds = driver.Create(path, self.x.size, self.y.size, 1, gdal.GDT_Float32)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(self.epsg)
    ds.SetProjection(srs.ExportToWkt())
    ds.SetGeoTransform(self.geoTransform)
    outband=ds.GetRasterBand(1)
    outband.SetStatistics(float(np.min(self.values)), float(np.max(self.values)),
                    float(np.average(self.values)), float(np.std(self.values)))
    outband.SetNoDataValue(-99999.0)
    if np.ma.is_masked(self.values):
        values =  np.ma.filled(self.values, -99999.0)
    else:
        values = self.values
    outband.WriteArray(values)


def apply_gaussian_filter(self, sigma, **kwargs):
    self.values = gaussian_filter(self.values, sigma, **kwargs)

def get_dx(self, meters=False):
    if meters is True and self.epsg==4326:
        return haversine((self.x[0], self.y[0]), (self.x[1], self.y[0]))*1000.
    else: return self.geoTransform[1] 
        
def get_dy(self, meters=False):
    if meters is True and self.epsg==4326:
        return haversine((self.x[0], self.y[1]), (self.x[0], self.y[0]))*1000.
    else: return -self.geoTransform[-1]
    
def get_extent(self, epsg=None):

    if epsg is None:
        epsg=self.epsg

    if self.epsg == epsg:
        return np.min(self.x), np.max(self.x), np.min(self.y), np.max(self.y)

    elif self.epsg==4326 and epsg != 4326:
        projection   = pyproj.Proj(init='epsg:'+str(epsg))
        min_x, min_y = projection(np.min(self.x), np.min(self.y))
        max_x, max_y = projection(np.max(self.x), np.max(self.y))
        return min_x, max_x, min_y, max_y
        
    elif self.epsg != 4326 and epsg == 4326:
        projection     = pyproj.Proj(init='epsg:'+str(self.epsg))
        minlon, minlat = projection(np.min(self.x), np.min(self.y), inverse=True)
        maxlon, maxlat = projection(np.max(self.x), np.max(self.y), inverse=True)
        return minlon, maxlon, minlat, maxlat

def resize_tile(self, dsfact):
    new_x = np.linspace(self.x[0], self.x[-1], self.x[::dsfact].size)
    new_y = np.linspace(self.y[0], self.y[-1], self.y[::dsfact].size)
    interp = RectBivariateSpline(self.x, self.y, self.values)
    new_z = interp(new_x, new_y)

    # set new pixel resolution in geoTransform
    res_x = (np.max(new_x)-np.min(new_x)) / (len(new_x)-1)
    res_y = ((np.min(new_y)-np.max(new_y)) / (len(new_y)-1))
   
    self.geoTransform = (self.geoTransform[0], res_x, self.geoTransform[2],
                            self.geoTransform[3], self.geoTransform[4], res_y)
    self.x = new_x
    self.y = new_y
    self.values = new_z

def get_xyz(self, epsg=None, include_invalid=False, path=None, radius=None, tansform=False):
    """
    Reshapes a DEM tile to a ndarray representing xyz coordinates.
    Output is a numpy array of shape (mx3) representing a "typical"
    'xyz' dataset or scattered point dataset.
    """
    x, y = np.meshgrid(self.x, self.y)
    x  = _x.reshape(x.size)
    y  = np.flipud(y.reshape(y.size))
    z = self.values.reshape(self.values.size)
    z = np.ma.filled(z, fill_value=np.nan)

    if epsg is None:
        epsg = self.epsg

    if self.epsg != epsg:
        target_proj = pyproj.Proj(init='epsg:{}'.format(epsg))
        if self.epsg!=4326:
            tile_proj = pyproj.Proj(init='epsg:{}'.format(self.epsg))
            _x, _y = pyproj.transform(tile_proj, target_proj, x, y)
        elif self.epsg==4326:
            _x, _y = target_proj(_x, _y)
        _x = np.asarray(_x).flatten()
        _y = np.asarray(_y).flatten()

    _xyz = np.vstack((_x, _y, z)).T
      
    if path is not None and include_invalid==False:
        idx, = np.where(np.logical_and(path.contains_points(_xyz[:,0:2], radius=radius), ~np.isnan(_xyz[:,2])))
    elif path is not None and include_invalid==True:
        idx, = np.where(path.contains_points(_xyz[:,0:2], radius=radius))
    elif path is None and include_invalid==False:
        idx, = np.where(~np.isnan(_xyz[:,2]))
    else:
        idx = np.arange(_xyz.shape[0])
    
    if transform==True:
        return _xyz[idx]  
    else:
        return np.vstack((x[idx], y[idx], z[idx])).T

def transform_to_epsg(self, epsg):
    self_proj = pyproj.Proj(init='epsg:{}'.format(self.epsg))
    target_proj = pyproj.Proj(init='epsg:{}'.format(epsg))
    x, y = pyproj.transform(self_proj, target_proj, self.x, self.y)
    self.x = np.asarray(x).flatten()
    self.y = np.asarray(y).flatten()
    self.geoTransform = ((self.x[0], (self.x[-1]-self.x[0])/(len(self.x)-1), 0, np.max(self.y), (np.min(self.y)-np.max(self.y))/(len(self.y)-1)))
    self.epsg = epsg

def concatenate_tiles(rootdir, extent, epsg, file_format):
    tile_list = list()
    for root, dirs, files in os.walk(rootdir):
        for file in files:
            if file.endswith(file_format):
                tile_list.append(root+"/"+file)
    xyz=list()
    for file in tile_list:
        # reading headers only for faster performance.
        geo  = gdal.Open(file)
        inproj = osr.SpatialReference()
        tile_epsg = inproj.GetAuthorityCode('PROJCS')
        if tile_epsg is None:
            tile_epsg = 4326
        geoTransform  = geo.GetGeoTransform()
        tile_min_x = geoTransform[0]
        tile_max_x = geoTransform[0] + geo.RasterXSize * geoTransform[1]
        tile_min_y = geoTransform[3] + geo.RasterYSize * geoTransform[5]
        tile_max_y = geoTransform[3]
   
        if tile_epsg != epsg:
            tile_proj = pyproj.Proj(init='epsg:{}'.format(tile_epsg))
            target_proj = pyproj.Proj(init='epsg:{}'.format(epsg))
            tile_min_x, tile_min_y = pyproj.transform(tile_proj, target_proj, tile_min_x, tile_min_y)
            tile_max_x, tile_max_y = pyproj.transform(tile_proj, target_proj, tile_max_x, tile_max_y)
        if (extent[0] <= tile_max_x and extent[1] >= tile_min_x) and \
            (extent[2] <= tile_max_y and extent[3] >= tile_min_y):
            tile = demtools.read_tile(file)
            xyz.append(tile.get_xyz(epsg=epsg))

    if len(xyz) > 0:
        return np.concatenate(tuple(xyz), axis=0)

def get_xyz_from_Path_instance(rootdir, Path_instance, epsg, file_format, radius=None, transform=False):

    tile_list = list()
    for root, dirs, files in os.walk(rootdir):
        for file in files:
            if file.endswith(file_format):# or file.endswith('_dem.img.xml'):
                tile_list.append(root+"/"+file)
    xyz = list()
    for file in tile_list:
        tile = demtools.read_tile(file)
        tile_path = tile.get_bbox_as_Path(epsg=epsg)
        if Path_instance.intersects_path(tile_path):
            xyz.append(tile.get_xyz(epsg=epsg, path=Path_instance, radius=radius, transform=transform))
    if len(xyz) > 0:
        return np.concatenate(tuple(xyz), axis=0)


def get_bbox_as_Path(self, **kwargs):
    target_epsg = kwargs.pop("epsg", self.epsg)
    if self.epsg != target_epsg:
        tile_proj = pyproj.Proj(init='epsg:{}'.format(self.epsg))
        target_proj = pyproj.Proj(init='epsg:{}'.format(target_epsg))
        tile_min_x, tile_min_y = pyproj.transform(tile_proj, target_proj, np.min(self.x), np.min(self.y))
        tile_max_x, tile_max_y = pyproj.transform(tile_proj, target_proj, np.max(self.x), np.max(self.y))
    else:
        tile_min_x = np.min(self.x)
        tile_min_y = np.min(self.y)
        tile_max_x = np.max(self.x)
        tile_max_y = np.max(self.y)
    path =  matplotlib.path.Path([(tile_min_x, tile_min_y),
                                            (tile_max_x, tile_min_y),
                                            (tile_max_x, tile_max_y),
                                            (tile_min_x, tile_max_y),
                                            (tile_min_x, tile_min_y)], closed=True)
    return path
        
    
    
        
        
#### Experimental functions begin here.

def scatter_to_geoTransform(xyz, geoTransform, epsg, shape, **kwargs):
    xpixels, ypixels = shape

    x = np.linspace(geoTransform[0], geoTransform[0] + xpixels*geoTransform[1], xpixels)
    y = np.linspace(geoTransform[3] + ypixels*geoTransform[5], geoTransform[3], ypixels)

    xt, yt = np.meshgrid(x, y)
    
    xt = xt.reshape(xt.size)
    yt = np.flipud(yt.reshape(yt.size))

    xyt = np.vstack((xt,yt)).T
    zt = griddata((xyz[:,0], xyz[:,1]), xyz[:,2], (xt, yt), method='linear', fill_value=np.nan)
    zt = zt.reshape(shape)
    return DEM.DEM(x, y, zt, geoTransform, epsg)


def build_DEM_from_files(rootdir, file_extension, extent, dx, dy):
    xyz = demtools.demutils.concatenate_tiles(rootdir, file_extension, extent)
    return demtools.demutils.scatter_to_tile(xyz, extent, dx, dy)


        

def deploy_to_PostGIS(self, dbname, **kwargs):
    user        = kwargs.pop('user', 'postgres')
    password    = kwargs.pop('password', True)
    host        = kwargs.pop('host', 'localhost')
    port        = kwargs.pop('port', 5432)
    overwrite   = kwargs.pop('overwrite', False)
    schema      = kwargs.pop('schema', 'raster')
    table       = kwargs.pop('table', 'raster')
    
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
        cur.execute('CREATE TABLE IF NOT EXISTS {}.{} (id SERIAL PRIMARY KEY, rast raster);'.format(schema, table))
        
        print("preparing array...")
        # array="["
        # for row in self.values:
            # array+="["
            # for number in row:
                # array+="{},"
            # array=array[-1]+"],"
        # array=array[-1]+"]"
        # array=np.ma.filled(self.values)
        # print(list(array))
       
        sql = "INSERT INTO {}.{} VALUES (".format(schema, table)
        sql+= "ST_SetValue("
        sql+= "ST_AddBand("
        sql+= "ST_MakeEmptyRaster({}, {}, {}, {}, {}, {}, {}, {}, srid={}),".format( \
              self.x.shape[0], self.y.shape[0], self.geoTransform[0], self.geoTransform[3],\
              self.geoTransform[1], self.geoTransform[5], self.geoTransform[2], self.geoTransform[4], self.epsg)
        sql+= "'64BF', NULL, NULL),"
        sql+= "1,1,1,ARRAY{}::double precision[][]));".format(str(np.ma.filled(self.values).tolist()))
        
        print("executing sql command...")
        cur.execute(sql)
        con.commit()
        con.close()
        
        
    
    
            # self.geoTransform[0] /* top left x */
            # self.geoTransform[1] /* w-e pixel resolution */
            # self.geoTransform[2] /* 0 */
            # self.geoTransform[3] /* top left y */
            # self.geoTransform[4] /* 0 */
            # self.geoTransform[5] /* n-s pixel resolution (negative value) */

# def _scatter_to_tile_finite_volume(xyz, extent, dx, dy):
    # """
    # Test implementation of Finite Volume interpolation for structured grids.
    # """
    # origin_x_m, origin_y_m, zone, letter = utm.from_latlon(extent[2], extent[0])
    # lat_step, lon_step = utm.to_latlon(origin_x_m + dx, origin_y_m + dy, zone, letter)
    # dx = lon_step - extent[0]
    # dy = lat_step - extent[2]
    # x = np.arange(extent[0], extent[1] + dx, step=dx)
    # y = np.arange(extent[2], extent[3] + dy, step=dy)
    # x, y = np.meshgrid(x,y)
    # m, n = x.shape
    # x = x.reshape(x.size)
    # y = np.flipud(y.reshape(y.size))
    # z = np.zeros(x.shape)
    # if progressbar:
        # widgets = ['[',Percentage(),']: ', Bar(marker='#',left='[',right=']'),
           # ' ', ETA()]
        # pbar = ProgressBar(widgets=widgets, maxval=x.size)
        # pbar.start()
    # for i in range(x.size):
        # idx, = np.where(np.logical_and(
            # np.logical_and(xyz[:,0]>=(xyz[i,0]-dx),xyz[:,0]<=(xyz[i,0]+dx)),
            # np.logical_and(xyz[:,1]>=(xyz[i,1]-dy),xyz[:,1]<=(xyz[i,1]+dy))))
        # z[i] = griddata((xyz[idx,0],xyz[idx,1]), xyz[idx,2], (x[i], y[i]), fill_value=-99999.0)
        # if progressbar:
            # pbar.update(i)
    # if progressbar:
        # pbar.finish()

    # x = x.reshape([m,n])
    # y = y.reshape([m,n])
    # z = z.reshape([m,n])
    # values = np.ma.masked_equal(z, -99999.0)
    # geoTransform = (extent[0], dx, 0, extent[-1], 0, -dy)
    # epsg = 4326
    # return x, y, values, geoTransform, epsg
    
# def scatter_to_tile(xyz, extent, dx, dy, epsg):

    # # We need to implement kernel interpolation. 
    # origin_x_m, origin_y_m, zone, letter = utm.from_latlon(extent[2], extent[0])
    # lat_step, lon_step = utm.to_latlon(origin_x_m + dx, origin_y_m + dy, zone, letter)
    # dx = lon_step - extent[0]
    # dy = lat_step - extent[2]
    # x = np.arange(extent[0], extent[1] + dx, step=dx)
    # y = np.arange(extent[2], extent[3] + dy, step=dy)
    # _x, _y = np.meshgrid(x, y)
    # m, n = _x.shape
    # _x = _x.reshape(_x.size)
    # _y = np.flipud(_y.reshape(_y.size))
    # z = griddata((xyz[:,0],xyz[:,1]), xyz[:,2], (_x,_y), fill_value=-99999.0)
    # z = z.reshape([m,n])
    # values = np.ma.masked_equal(z, -99999.0)
    # geoTransform = (extent[0], dx, 0, extent[-1], 0, -dy)
    # epsg = 4326
    # return DEM.DEM(x, y, values, geoTransform, epsg)