from __future__ import absolute_import, division, print_function
import numpy as np
from . import demutils
from . import demplots


def read_tile(path):
    """
    Reads DEM files of any format supported by GDAL.
        
    args:
        path (string): Path to DEM file.
    
    """
    return demutils.read_tile(path)


def concatenate_tiles(rootdir, extent, epsg, file_format):
    """
    Searches a directory recursively for all the DEM files that contain data within
    the bounding box defined by extent and returns a [D,3] numpy array of the concatenated
    data. Useful for passing data to interpolators and creating seamless DEM's from multiple
    data sources.
    
    args:
        rootdir (string)        : Top level directory from which to initilize the search.
        file_extension (string) : Will only read metadata from file ending with this extension.
        extent (list)           : Bounding box of the form [min_x, max_x, min_y, max_y].
        epsg                    : EPSG code for the coordinates given in extent.
    returns:
        Numpy array of dimension [D,3] containing the concatenation of all the DEM data for under 
        specified rootdir, bounded by specified extent and in epsg coordinates.
    """
    return demutils.concatenate_tiles(rootdir, extent, epsg, file_format)

def scatter_to_geoTransform(xyz, geoTransform, epsg, shape, **kwargs):
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
    return demutils.scatter_to_geoTransform(xyz, geoTransform, epsg, shape, **kwargs)

    
def build_DEM_from_files(rootdir, extent, dx, dy, file_format):
    """
    ** PROTOTYPE FUNCTION ** 
    ** WILL BE MODIFIED IN NON_BACKWARDS COMPATIBLE WAY **

    Searches recursively a directory for all DEM files that contain data within a given 
    bounding box, and generates a seamless DEM object of resolution dx, dy.
    
    args:
        rootdir (string)        : Top level directory from which to initilize the search.
        file_extension (string) : Will only read metadata from file ending with this extension.
        extent (list)           : Bounding box of the form [min_lon, max_lon, min_lat, max_lat].
                                  Requires lonlat coordinates regardless of the DEM tile projection.
        dx                      : Output resolution in x direction given in meters.
        dy                      : Output resolution in y direction given in meters.
    
    Notes:
        This function will be extended in the future to accept inputs based on epsg code
        without the hardcoded latlon system it currently implements.
    """
    return demutils.build_DEM_from_files(rootdir, extent, dx, dy, file_format)

def get_xyz_from_Path_instance(rootdir, Path, epsg, file_format, radius=None, transform=False):
    """
    Will return the data on the same EPSG code as the EPSG code of the Path instance specified.
    """
    return demutils.get_xyz_from_Path_instance(rootdir, Path, epsg, file_format, radius, transform)


class DEM(object):
    """
    Main DEM class. Contains DEM data attributes and multiple functions to manipulate
    and export DEM data.

    args:
    -----
        x (numpy array) 
    """
    def __init__(self, x, y, values, geoTransform, epsg, vertical_datum=None):
        """
        Initiliazes the DEM object from arrays.
        
        args:
            x (numpy array, flattened) : 1D numpy array of x-coordinate values.
            
        
        """
        self.x            = x
        self.y            = y
        self.values       = values
        self.epsg         = epsg
        self.geoTransform = geoTransform
        self.vertical_datum = vertical_datum
        
    def __sub__(self, other):
        values = self.values - other.values
        return DEMdiff(self.y, self.y, values, self.geoTransform, self.epsg)
               
    def make_plot(self, axes=None, vmin=None, vmax=None, title=None, colors=256):
        return demplots.plot_bathy(self, axes, vmin, vmax, title, colors)
    
    def get_extent(self, epsg=None):
        return demutils.get_extent(self, epsg)
        
    def get_dx(self, meters=False):
        return demutils.get_dx(self, meters)
        
    def get_dy(self, meters=False):
        return demutils.get_dy(self, meters)

    def export_to_file(self, path, driver='Gtiff'):
        """
        Exports DEM object to raster file using GDAL.
        Uses GeoTiff as default. 

        Args:
        -----

            path (str) : Path to export file to.

        kwargs:
        -------

            driver (str) <'Gtiff'> : Specifies export file format supported by GDAL.

        """
        demutils.export_to_file(self, path, driver)
    
    def get_resolution(self):
        return demutils.resolution_in_meters(self)

    def get_xyz(self, **kwargs):
        """
        Returns a [D,3] numpy array with columns x, y, values.
        kwargs:
            epsg
        """
        return demutils.get_xyz(self, **kwargs)

    def downsample_by_factor(self, dsfact):
        demutils.resize_tile(self, dsfact)
        
    def apply_gaussian_filter(self, sigma, **kwargs):
        return demutils.apply_gaussian_filter(self, sigma, **kwargs)

    def apply_circular_filter(self, radius):
        demutils.circular_filter(self, radius)
        
    def apply_circular_pixel_filter(self, radius):
        demutils.circular_pixel_filter(self, radius)
    
    def apply_rbf_filter(self):
        pass   

    def write_to_PostGIS(self, **kwargs):
        raise NotImplementedError("Coming soon!")

    def relative_RMSE(self, other):
        return np.mean(np.sqrt(np.square(self.values - other.values)))
        
    def deploy_to_PostGIS(self, **kwargs):
        demutils.deploy_to_PostGIS(self, **kwargs)
        
    def transform_to_epsg(self, epsg):
        demutils.transform_to_epsg(self, epsg)
    
    def get_bbox_as_Path(self, **kwargs):
        return demutils.get_bbox_as_Path(self, **kwargs)
        
class DEMdiff(DEM):
    def __init__(self, x, y, diff, geoTransform, epsg):
        super(DEMdiff, self).__init__(x, y, diff, geoTransform, epsg)
    
    def make_plot(self, axes=None, vmin=None, vmax=None, title=None):
        return demplots.plot_diff(self, axes, vmin, vmax, title)

