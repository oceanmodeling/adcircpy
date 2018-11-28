import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pyproj
from osgeo import osr, gdal
from AdcircPy.core import FixPointNormalize
gdal.UseExceptions()

class DEM(object):
  """
  Main DEM class. Contains DEM data attributes and multiple functions to manipulate
  and export DEM data.

  args:
  -----
      x (numpy array) 
  """
  def __init__(self, x, y, values, geoTransform, epsg, vertical_datum=None):
    self.x            = x
    self.y            = y
    self.values       = values
    self.epsg         = epsg
    self.geoTransform = geoTransform
    self.vertical_datum = vertical_datum

  def __sub__(self, other):
    values = self.values - other.values
    return DEMdiff(self.y, self.y, values, self.geoTransform, self.epsg)

  @classmethod
  def read_tile (cls, path, epsg=None, vertical_datum=None):
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
    if vertical_datum is None:
      vertical_datum = inproj.GetAttrValue('datum')
    if epsg is None:
      if inproj.IsGeographic() == 1:
        epsg = 4326
      else:
        epsg = inproj.GetAttrValue('PROJCS|AUTHORITY', 1)
        # Hack to force identification of WGS84 on certain files.
        if epsg is None and geo.GetGeoTransform()[2]<10**-4:
          epsg=4326
        elif epsg is None:
          raise Exception("Could not auto identify tile EPSG. Please set EPSG using epsg kwarg. \nFore reference, WKT is {}".format(inproj))
    geoTransform  = geo.GetGeoTransform()
    x = np.linspace(geoTransform[0], 
                    geoTransform[0] + geo.RasterXSize * geoTransform[1],
                    num = geo.RasterXSize)
    y = np.linspace(geoTransform[3] + geo.RasterYSize * geoTransform[5],
                    geoTransform[3],
                    num = geo.RasterYSize)
    values = geo.ReadAsArray()
    values = np.ma.masked_equal(values, geo.GetRasterBand(1).GetNoDataValue())
    return cls(x, y, values, geoTransform, epsg, vertical_datum)
             
  @classmethod
  def interpolate_from_xyz(cls, xyz, geoTransform, epsg, shape, **kwargs):
    """
    Interpolates numpy array of the form [D,3] into bounding box specified by extent
    and resolution dx, dy.
    
    args:
        xyz (numpy array) : Dimensions are [D,3].
        geoTansform       : 
        epsg              : EPSG code of xyz data and geoTransform data (must match).
        shape             : Shape of resulting DEM.
    
    return:
        DEM object in latlon coordinates (epsg 4326).
        
    Notes:
        This function will be extended in the future to accept inputs based on epsg code
        without the hardcoded latlon system it currently implements.
    
    """
    xpixels, ypixels = shape

    x = np.linspace(geoTransform[0], geoTransform[0] + xpixels*geoTransform[1], xpixels)
    y = np.linspace(geoTransform[3] + ypixels*geoTransform[5], geoTransform[3], ypixels)

    xt, yt = np.meshgrid(x, y)
    
    xt = xt.reshape(xt.size)
    yt = np.flipud(yt.reshape(yt.size))

    xyt = np.vstack((xt,yt)).T
    zt = griddata((xyz[:,0], xyz[:,1]), xyz[:,2], (xt, yt), method='linear', fill_value=np.nan)
    zt = zt.reshape(shape)
    return cls(x, y, zt, geoTransform, epsg)
          
  @classmethod
  def get_xyz_from_Path_instance(cls, rootdir, Path_instance, epsg, transform=False):
    tile_list = list()
    for root, dirs, files in os.walk(rootdir):
      for file in files:
        tile_list.append(root+"/"+file)
    xyz = list()
    for file in tile_list:
      # try:
        tile = cls.read_tile(file)
        tile_path = tile.get_bbox_as_Path(epsg=epsg)
        if Path_instance.intersects_path(tile_path):
          xyz.append(tile.get_xyz(epsg=epsg, path=Path_instance, transform=transform))
      # except:
      #   continue
    if len(xyz) > 0:
      return np.concatenate(tuple(xyz), axis=0)

  def make_plot(self, axes=None, vmin=None, vmax=None, title=None, colors=256, show=False):
      if axes is None:
          fig = plt.figure()
          axes = fig.add_subplot(111)

      if vmin is None:
          vmin = np.ceil(np.min(self.values))
      if vmax is None:
          vmax = np.floor(np.max(self.values))
      
      min_z = vmin
      max_z = vmax        

      if max_z < 0.:
          mlevel = np.mean(self.values)
          cmap = plt.cm.seismic
          col_val = 0.5
          levels=np.linspace(min_z, max_z, colors)        
      else:
          wet_count = int(np.floor(float(colors) * (float((self.values < 0.).sum())/float(self.values.size))))
          col_val = float(wet_count)/float(colors)
          dry_count = colors - wet_count
          colors_undersea = plt.cm.bwr(np.linspace(1., 0., wet_count))
          colors_land = plt.cm.terrain(np.linspace(0.25, 1., dry_count))
          colors = np.vstack((colors_undersea, colors_land))
          cmap = LinearSegmentedColormap.from_list('cut_terrain', colors)
          wlevels = np.linspace(min_z, 0.0, wet_count, endpoint=False)
          dlevels = np.linspace(0.0, max_z, dry_count)
          levels = np.hstack((wlevels, dlevels))
      norm = FixPointNormalize(sealevel=0.0, vmax=max_z, vmin=min_z, col_val=col_val)
      ax = axes.pcolormesh(self.x, self.y, np.flipud(self.values), cmap=cmap, norm=norm)
      axes.axis('scaled')
      if title is not None:
          axes.set_title(title)
      mappable = ScalarMappable(cmap=cmap)
      mappable.set_array([])
      mappable.set_clim(min_z, max_z)
      divider = make_axes_locatable(axes)
      cax = divider.append_axes("bottom", size="2%", pad=0.5)
      cbar = plt.colorbar(mappable, cax=cax,
                          label=r'elevation [m]',
                          extend='both', orientation='horizontal')

      cbar.set_ticks([vmin, vmin + col_val *(vmax-vmin), vmax])
      cbar.set_ticklabels([np.around(vmin, 2), 0.0, np.around(vmax, 2)])
      if show==True:
        plt.show()
      
      return axes
  
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
      
  def get_dx(self, meters=False):
    if meters is True and self.epsg==4326:
      return haversine((self.x[0], self.y[0]), (self.x[1], self.y[0]))*1000.
    else:
      return self.geoTransform[1] 

  def get_dy(self, meters=False):
    if meters is True and self.epsg==4326:
      return haversine((self.x[0], self.y[1]), (self.x[0], self.y[0]))*1000.
    else:
      return -self.geoTransform[-1]
  
  def dump(self, path, driver='Gtiff'):
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
  
  def get_xyz(self, epsg=None, path=None, transform=None, include_invalid=False):
    """
    Reshapes a DEM tile to a ndarray representing xyz coordinates.
    Output is a numpy array of shape (mx3) representing a "typical"
    'xyz' dataset or scattered point dataset.
    """
    x, y = np.meshgrid(self.x, self.y)
    x  = x.reshape(x.size)
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
        _x, _y = target_proj(x, y)
      _x = np.asarray(_x).flatten()
      _y = np.asarray(_y).flatten()
      _xyz = np.vstack((_x, _y, z)).T
    elif self.epsg==epsg:
      _xyz = np.vstack((x,y,z)).T
    if path is not None and include_invalid==False:
      idx, = np.where(np.logical_and(path.contains_points(_xyz[:,0:2]), ~np.isnan(_xyz[:,2])))
    elif path is not None and include_invalid==True:
      idx, = np.where(path.contains_points(_xyz[:,0:2]))
    elif path is None and include_invalid==False:
      idx, = np.where(~np.isnan(_xyz[:,2]))
    elif path is None and include_invalid==True:
      idx = np.arange(_xyz.shape[0])
    else:
      raise Exception("This exception should be unreachable.")
    if transform==True:
      return _xyz[idx,:]
    else:
      return np.vstack((x[idx], y[idx], z[idx])).T

  def downsample_by_factor(self, dsfact):
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
      
  def apply_gaussian_filter(self, sigma, **kwargs):
    """ """
    self.values = gaussian_filter(self.values, sigma, **kwargs)

  def apply_circular_filter(self, radius):
    """
    Input radius in meters.
    """
    yres=self.get_dy(meters=True)
    xres=self.get_dx(meters=True)
    num = 0.
    while num < radius:
      num = num + xres
    x_radius = num - xres
    num = 0.
    while num < radius:
      num = num + yres
    y_radius = num - yres
    x = np.arange(-x_radius, x_radius+xres, xres)
    y = np.arange(-y_radius, y_radius+yres, yres)
    y,x = np.ogrid[-radius:radius+yres:yres, -radius:radius+xres:xres]
    mask = x**2 + y**2 <= radius**2
    kernel = np.zeros((mask.shape[0], mask.shape[1]))
    kernel[mask] = 1
    self.values = generic_filter(self.values, np.mean, footprint=kernel)
 
  def apply_circular_pixel_filter(self, radius):
    """
    Input radius in meters.
    """
    ncells = int(radius // resolution_in_meters(self))
    kernel = np.zeros((2*ncells+1, 2*ncells+1))
    y,x = np.ogrid[-ncells:ncells+1, -ncells:ncells+1]
    mask = x**2 + y**2 <= ncells**2
    kernel[mask] = 1
    self.values = generic_filter(self.values, np.mean, footprint=kernel)
  
  def relative_RMSE(self, other):
    """ """
    return np.mean(np.sqrt(np.square(self.values - other.values)))
      
  def transform_to_epsg(self, epsg):
    self_proj = pyproj.Proj(init='epsg:{}'.format(self.epsg))
    target_proj = pyproj.Proj(init='epsg:{}'.format(epsg))
    x, y = pyproj.transform(self_proj, target_proj, self.x, self.y)
    self.x = np.asarray(x).flatten()
    self.y = np.asarray(y).flatten()
    self.geoTransform = ((self.x[0], (self.x[-1]-self.x[0])/(len(self.x)-1), 0, np.max(self.y), (np.min(self.y)-np.max(self.y))/(len(self.y)-1)))
    self.epsg = epsg
  
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
    path = Path([(tile_min_x, tile_min_y),
                                            (tile_max_x, tile_min_y),
                                            (tile_max_x, tile_max_y),
                                            (tile_min_x, tile_max_y),
                                            (tile_min_x, tile_min_y)], closed=True)
    return path

  def deploy_to_PostGIS(self, dbname, **kwargs):
    raise NotImplementedError
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
      
  def __build_DEM_from_files(rootdir, extent, dx, dy, file_format):
    """ placeholder, might be removed """
    raise NotImplementedError


class DEMdiff(DEM):
  def __init__(self, x, y, diff, geoTransform, epsg):
    super(DEMdiff, self).__init__(x, y, diff, geoTransform, epsg)
  
  def make_plot(self, axes=None, vmin=None, vmax=None, title=None):
    if axes is None:
      fig = plt.figure()
      axes = fig.add_subplot(111)
    if vmin is None:
      vmin = np.ceil(np.min(self.values))
    if vmax is None:
      vmax = np.floor(np.max(self.values))
    cmap='seismic'
    norm = FixPointNormalize(sealevel=0, vmax=vmax, vmin=vmin, col_val=0.5)
    axes.imshow(self.values, origin='upper',
                                            extent=self.get_extent(),
                                            cmap=cmap,
                                            norm=norm)

    if title is not None:
      axes.set_title(title)

    mappable = matplotlib.cm.ScalarMappable(cmap=cmap)
    mappable.set_array([])
    mappable.set_clim(vmin, vmax)
    divider = make_axes_locatable(axes)
    cax = divider.append_axes("bottom", size="2%", pad=0.5)
    cbar = plt.colorbar(mappable, cax=cax,
                        label=r'elevation [$\Delta$m]',
                        extend='both', orientation='horizontal')
    cbar.set_ticks([vmin, vmin + 0.5*(vmax-vmin), vmax])
    cbar.set_ticklabels([np.around(vmin, 2), 0.0, np.around(vmax, 2)])
    return axes
        

