import numpy as np
from matplotlib.path import Path
from matplotlib.transforms import Bbox
from osgeo import osr, gdal


class gdal_tools(object):
    """
    From gdal docs:
        adfGeoTransform[0] /* top left x */
        adfGeoTransform[1] /* w-e pixel resolution */
        adfGeoTransform[2] /* 0 */
        adfGeoTransform[3] /* top left y */
        adfGeoTransform[4] /* 0 */
        adfGeoTransform[5] /* n-s pixel resolution (negative value) */
    """

    @classmethod
    def Warp(cls, Dataset, filename='', format='VRT', **kwargs):
        return gdal.Warp(filename, Dataset, format=format, **kwargs)

    @classmethod
    def get_Bbox(cls, Dataset, SpatialReference=None, Path=False):
        gt = Dataset.GetGeoTransform()
        cols = Dataset.RasterXSize
        rows = Dataset.RasterYSize
        xarr = [0, cols]
        yarr = [0, rows]
        ext = list()  # There should be a more succint method.
        for px in xarr:
            for py in yarr:
                x = gt[0]+(px*gt[1])+(py*gt[2])
                y = gt[3]+(px*gt[4])+(py*gt[5])
                ext.append([x, y])
            yarr.reverse()
        ext = np.asarray(ext)
        xmin = np.min(ext[:, 0])
        xmax = np.max(ext[:, 0])
        ymin = np.min(ext[:, 1])
        ymax = np.max(ext[:, 1])
        bbox = Bbox([[xmin, ymin], [xmax, ymax]])
        if SpatialReference is not None:
            if isinstance(SpatialReference, int):
                epsg = SpatialReference
                SpatialReference = osr.SpatialReference()
                SpatialReference.ImportFromEPSG(epsg)
            else:
                assert isinstance(SpatialReference, osr.SpatialReference)
            inSpatialRef = cls.get_SpatialReference(Dataset)
            CoordinateTransform = osr.CoordinateTransformation(
                                                            inSpatialRef,
                                                            SpatialReference)
            p = CoordinateTransform.TransformPoints([[bbox.xmin, bbox.ymin],
                                                    [bbox.xmax, bbox.ymax]])
            bbox = Bbox([[p[0][0], p[0][1]], [p[1][0], p[1][1]]])
        if Path:
            bbox = cls.bbox_to_Path(bbox)
        return bbox

    @classmethod
    def get_xyz(cls, Dataset, SpatialReference=None, bbox=None):
        x, y, z = cls.get_arrays(Dataset, SpatialReference)
        x, y = np.meshgrid(x, y)
        x = x.flatten()
        y = y.flatten()
        z = z.flatten()
        xyz = np.vstack([x, y, z]).T
        if bbox is not None:
            idxs, = np.where(np.logical_and(
                             np.logical_and(xyz[:, 0] >= bbox.xmin,
                                            xyz[:, 0] <= bbox.xmax),
                             np.logical_and(xyz[:, 1] >= bbox.ymin,
                                            xyz[:, 1] <= bbox.ymax)))
            xyz = xyz[idxs, :]
        return xyz

    @classmethod
    def get_arrays(cls, Dataset, SpatialReference=None):
        xmin, dx, _, ymax, _, dy = cls.get_GeoTransform(Dataset,
                                                        SpatialReference)
        ymin = ymax + Dataset.RasterYSize*dy
        xmax = xmin + Dataset.RasterXSize*dx
        x = np.linspace(xmin, xmax, Dataset.RasterXSize)
        y = np.linspace(ymin, ymax, Dataset.RasterYSize)
        z = np.flipud(Dataset.ReadAsArray())
        return x, y, z

    @classmethod
    def get_resolution(cls, Dataset, SpatialReference=None):
        _, dx, _, _, _, dy = cls.get_GeoTransform(Dataset, SpatialReference)
        return dx, dy

    @classmethod
    def get_GeoTransform(cls, Dataset, SpatialReference=None):
        xmin, xres, A, ymax, B, yres = Dataset.GetGeoTransform()
        if SpatialReference is not None:
            bbox = cls.get_Bbox(Dataset, SpatialReference)
            xmin = bbox.xmin
            xres = (bbox.xmax-xmin)/Dataset.RasterXSize
            ymax = bbox.ymax
            yres = (bbox.ymin-ymax)/Dataset.RasterYSize
        return xmin, xres, A, ymax, B, yres

    @staticmethod
    def bbox_to_Path(bbox):
        bbox = Bbox(bbox)
        return Path([[bbox.xmin, bbox.ymin],
                     [bbox.xmax, bbox.ymin],
                     [bbox.xmax, bbox.ymax],
                     [bbox.xmin, bbox.ymax],
                     [bbox.xmin, bbox.ymin]])

    @staticmethod
    def get_SpatialReference(Dataset):
        return osr.SpatialReference(wkt=Dataset.GetProjection())


resample = gdal_tools.Warp
get_Bbox = gdal_tools.get_Bbox
get_xyz = gdal_tools.get_xyz
get_arrays = gdal_tools.get_arrays
get_resolution = gdal_tools.get_resolution
get_GeoTransform = gdal_tools.get_GeoTransform
bbox_to_Path = gdal_tools.bbox_to_Path
get_SpatialReference = gdal_tools.get_SpatialReference
