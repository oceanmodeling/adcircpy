# gobal imports
import numpy as np
import pyproj
from matplotlib.path import Path
from scipy.interpolate import RectBivariateSpline
from haversine import haversine
from osgeo import gdal, osr


class _BaseRaster(object):

    def __init__(self, x, y, values, GeoTransform, SpatialReference):
        self._x = x
        self._y = y
        self._values = values
        self._SpatialReference = SpatialReference
        self._GeoTransform = GeoTransform
        super(_BaseRaster, self).__init__()

    def __sub__(self, other):
        # importing here to avoid cyclical dependency
        from AdcircPy.DEM._RasterDiff import _RasterDiff
        return _RasterDiff(self.x, self.y, self.values - other.values,
                           self.GeoTransform, self.SpatialReference)

    def get_extent(self, epsg=None):
        if epsg is None:
            epsg = self.epsg
        if self.epsg == epsg:
            return np.min(self.x), np.max(self.x), np.min(self.y),
            np.max(self.y)
        elif self.epsg == 4326 and epsg != 4326:
            projection = pyproj.Proj(init='epsg:'+str(epsg))
            min_x, min_y = projection(np.min(self.x), np.min(self.y))
            max_x, max_y = projection(np.max(self.x), np.max(self.y))
            return min_x, max_x, min_y, max_y
        elif self.epsg != 4326 and epsg == 4326:
            projection = pyproj.Proj(init='epsg:'+str(self.epsg))
            minlon, minlat = projection(np.min(self.x), np.min(self.y),
                                        inverse=True)
            maxlon, maxlat = projection(np.max(self.x), np.max(self.y),
                                        inverse=True)
            return minlon, maxlon, minlat, maxlat

    def get_dx(self, meters=False):
        if meters is True and self.epsg == 4326:
            return haversine((self.x[0], self.y[0]),
                             (self.x[1], self.y[0]))*1000.
        else:
            return self.geoTransform[1] 

    def get_dy(self, meters=False):
        if meters is True and self.epsg == 4326:
            return haversine((self.x[0], self.y[1]),
                             (self.x[0], self.y[0]))*1000.
        else:
            return -self.geoTransform[-1]

    def dump(self, path, driver='Gtiff'):
        driver = gdal.GetDriverByName(driver)
        ds = driver.Create(path, self.values.shape[1], self.values.shape[0], 1,
                           gdal.GDT_Float32)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(self.epsg)
        ds.SetProjection(srs.ExportToWkt())
        ds.SetGeoTransform(self.geoTransform)
        outband = ds.GetRasterBand(1)
        outband.SetStatistics(float(np.min(self.values)),
                              float(np.max(self.values)),
                              float(np.average(self.values)),
                              float(np.std(self.values)))
        outband.SetNoDataValue(-99999.0)
        if np.ma.is_masked(self.values):
            values = np.ma.filled(self.values, -99999.0)
        else:
            values = self.values
        outband.WriteArray(values)

    def get_xyz(self, epsg=None, path=None, transform=None,
                include_invalid=False):
        """
        Reshapes a DEM tile to a ndarray representing xyz coordinates.
        Output is a numpy array of shape (mx3) representing a "typical"
        'xyz' dataset or scattered point dataset.
        """
        x, y = np.meshgrid(self.x, self.y)
        x = x.reshape(x.size)
        y = np.flipud(y.reshape(y.size))
        z = self.values.reshape(self.values.size)
        z = np.ma.filled(z, fill_value=np.nan)
        if epsg is None:
            epsg = self.epsg
        if self.epsg != epsg:
            target_proj = pyproj.Proj(init='epsg:{}'.format(epsg))
            if self.epsg != 4326:
                tile_proj = pyproj.Proj(init='epsg:{}'.format(self.epsg))
                _x, _y = pyproj.transform(tile_proj, target_proj, x, y)
            elif self.epsg == 4326:
                _x, _y = target_proj(x, y)
            _x = np.asarray(_x).flatten()
            _y = np.asarray(_y).flatten()
            _xyz = np.vstack((_x, _y, z)).T
        elif self.epsg == epsg:
            _xyz = np.vstack((x, y, z)).T
        if path is not None and include_invalid is False:
            idx, = np.where(np.logical_and(path.contains_points(_xyz[:, 0:2]),
                                           ~np.isnan(_xyz[:, 2])))
        elif path is not None and include_invalid is True:
            idx, = np.where(path.contains_points(_xyz[:, 0:2]))
        elif path is None and include_invalid is False:
            idx, = np.where(~np.isnan(_xyz[:, 2]))
        elif path is None and include_invalid is True:
            idx = np.arange(_xyz.shape[0])
        else:
            raise Exception("This exception should be unreachable.")
        if transform is True:
            return _xyz[idx, :]
        else:
            return np.vstack((x[idx], y[idx], z[idx])).T

    def downsample_by_factor(self, dsfact):
        """
        Downsamples DEM by removing each dsfact points, effectively
        lowering resolution.
        """
        new_x = np.linspace(self.x[0], self.x[-1], self.x[::dsfact].size)
        new_y = np.linspace(self.y[0], self.y[-1], self.y[::dsfact].size)
        interp = RectBivariateSpline(self.x, self.y, self.values)
        new_z = interp(new_x, new_y)
        # set new pixel resolution in geoTransform
        res_x = (np.max(new_x)-np.min(new_x)) / (len(new_x)-1)
        res_y = ((np.min(new_y)-np.max(new_y)) / (len(new_y)-1))
        self._geoTransform = (self.geoTransform[0], res_x,
                              self.geoTransform[2], self.geoTransform[3],
                              self.geoTransform[4], res_y)
        self._x = new_x
        self._y = new_y
        self._values = new_z

    def transform_to_epsg(self, epsg):
        self_proj = pyproj.Proj(init='epsg:{}'.format(self.epsg))
        target_proj = pyproj.Proj(init='epsg:{}'.format(epsg))
        x, y = pyproj.transform(self_proj, target_proj, self.x, self.y)
        self._x = np.asarray(x).flatten()
        self._y = np.asarray(y).flatten()
        self._geoTransform = ((self.x[0], (self.x[-1]-self.x[0]) /
                               (len(self.x)-1), 0, np.max(self.y),
                               (np.min(self.y)-np.max(self.y)) /
                               (len(self.y)-1)))
        self._epsg = epsg

    def get_bbox_as_Path(self, **kwargs):
        target_epsg = kwargs.pop("epsg", self.epsg)
        if self.epsg != target_epsg:
            tile_proj = pyproj.Proj(init='epsg:{}'.format(self.epsg))
            target_proj = pyproj.Proj(init='epsg:{}'.format(target_epsg))
            tile_min_x, tile_min_y = pyproj.transform(tile_proj, target_proj,
                                                      np.min(self.x),
                                                      np.min(self.y))
            tile_max_x, tile_max_y = pyproj.transform(tile_proj, target_proj,
                                                      np.max(self.x),
                                                      np.max(self.y))
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

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def values(self):
        return self._values

    @property
    def GeoTransform(self):
        return self._GeoTransform

    @property
    def SpatialReference(self):
        return self._SpatialReference

    @property
    def xyz(self):
        x, y = np.meshgrid(self.x, self.y)
        x = x.flatten()
        y = y.flatten()
        z = self.values.flatten()
        return np.vstack([x, y, z]).T
