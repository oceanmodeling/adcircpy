import os
import numpy as np
from osgeo import osr, gdal
from AdcircPy.DEM import _BaseDEM, _DEMDiff
from AdcircPy.core import _PostGIS

class DEM(_BaseDEM, _PostGIS):
  """
  Main DEM class. Contains DEM data attributes and multiple functions to manipulate
  and export DEM data.

  args:
  -----
      x (numpy array) 
  """
  def __init__(self, x, y, values, geoTransform, epsg, vertical_datum, Dataset=None):
    super(DEM, self).__init__(x, y, values, geoTransform, epsg)
    self._vertical_datum = vertical_datum
    self._Dataset = Dataset

  def __sub__(self, other):
    values = self.values - other.values
    return _DEMDiff(self.y, self.y, values, self.geoTransform, self.epsg)

  @property
  def vertical_datum(self):
    return self._vertical_datum

  @property
  def Dataset(self):
    return self._Dataset

  @classmethod
  def from_file(cls, path):
    return cls.__init_dataset(gdal.Open(path))

  @classmethod
  def from_PostGIS(cls, dbname, schema, table, user='postgres', password=None, host='localhost', port=5432, where=None, mode=None):
    if password is None:
      password = getpass.getpass(prompt='Enter Postgres password: ')
    conn_string = cls._get_PG_conn_string(dbname, host, port, user, password)
    dataset = gdal.Open(conn_string)
    subdatasets = dataset.GetSubDatasets()
    for dataset in subdatasets:
      if "schema='{}'".format(schema) in dataset[0] and "table='{}'".format(table) in dataset[0]:
        break
    Dataset = gdal.Open(dataset[0])
    return cls.__init_dataset(Dataset)
  
  @classmethod
  def from_xyz(cls, xyz, geoTransform, epsg, shape, method='linear', fill_value=np.nan, **kwargs):
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
    zt = griddata((xyz[:,0], xyz[:,1]), xyz[:,2], (xt, yt), method=method, fill_value=fill_value)
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

  def make_plot(self, axes=None, vmin='auto', vmax='auto', title=None, show=False, figsize=None):
    self._init_make_plot('topobathy', axes, vmin, vmax, figsize)
    self._ax = self._axes.pcolormesh(self.x, self.y, np.flipud(self.values), cmap=self._cmap, norm=self._norm)
    self._finalize_make_plot(title, show, 'neither', cbar_label=r'elevation [m]')

  def apply_gaussian_filter(self, sigma, **kwargs):
    """ """
    self._values = gaussian_filter(self.values, sigma, **kwargs)

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
      
  @classmethod
  def __init_dataset(cls, dataset):
    inproj = osr.SpatialReference(wkt=dataset.GetProjection())
    vertical_datum = inproj.GetAttrValue('datum')
    if inproj.IsGeographic() == 1:
      epsg = 4326
    else:
      epsg = inproj.GetAttrValue('PROJCS|AUTHORITY', 1)
      # Hack to force identification of WGS84 on certain files.
      if epsg is None and dataset.GetGeoTransform()[2]<10**-4:
        epsg=4326
    geoTransform  = dataset.GetGeoTransform()
    x = np.linspace(geoTransform[0], 
                    geoTransform[0] + dataset.RasterXSize * geoTransform[1],
                    num = dataset.RasterXSize)
    y = np.linspace(geoTransform[3] + dataset.RasterYSize * geoTransform[5],
                    geoTransform[3],
                    num = dataset.RasterYSize)
    values = dataset.ReadAsArray()
    values = np.ma.masked_equal(values, dataset.GetRasterBand(1).GetNoDataValue())
    return cls(x, y, values, geoTransform, epsg, vertical_datum, dataset) 

if __name__ == '__main__':
  tile = DEM.from_file(os.getenv('TESTDEM1'))
  tile.make_plot(show=True)

