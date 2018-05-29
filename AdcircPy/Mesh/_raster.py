import numpy as np
from matplotlib.path import Path
from scipy.interpolate import griddata
from AdcircPy import DEM


def rasterize_to_geoTransform(self, geoTransform, shape, **kwargs):
    """
    Converts an ADCIRC mesh object into a ADCIRC raster object.
    With this

    Parameters
    ----------
    geoTransform : tuple
        geoTransform[0] /* top left x */
        geoTransform[1] /* west-east pixel resolution */
        geoTransform[2] /* 0 */
        geoTransform[3] /* top left y */
        geoTransform[4] /* 0 */
        geoTransform[5] /* north-south pixel resolution (negative value) */

    xpixels : int
        Description of arg2

    ypixels : int

    Returns
    -------
    raster object
        Can be used to export to geoTif, and other operations such as filtering.
    """
    epsg = kwargs.pop("epsg", self.epsg)
    padding = kwargs.pop("padding", None)
    
    xpixels, ypixels = shape

    x = np.linspace(geoTransform[0], geoTransform[0] + xpixels*geoTransform[1], xpixels)
    y = np.linspace(geoTransform[3] + ypixels*geoTransform[5], geoTransform[3], ypixels)

    xt, yt = np.meshgrid(x, y)
    
    xt = xt.reshape(xt.size)
    yt = np.flipud(yt.reshape(yt.size))

    xyt = np.vstack((xt,yt)).T

    # create path object of target bounding box
    bbox_path = Path([(np.min(xt), np.min(yt)),
                    (np.max(xt), np.min(yt)),
                    (np.max(xt), np.max(yt)),
                    (np.min(xt), np.max(yt)),
                    (np.min(xt), np.min(yt))], closed=True)

    # interpolate mesh information to bounding box grid.
    xyz = self.get_xyz(extent=bbox_path)#, radius=0.2)
    zt = griddata((xyz[:,0], xyz[:,1]), xyz[:,2], (xt, yt), method='linear', fill_value=np.nan)
    
    # Generate boundary masks.
    outerBoundary = self.build_outer_polygon()
    outerBoundary = outerBoundary.clip_to_bbox((np.min(xt), np.min(yt), np.max(xt), np.max(yt)))
    mask = np.logical_or(np.isin(zt, np.nan), ~outerBoundary.contains_points(xyt))
    zt = np.ma.masked_array(zt, mask)
    
    innerBoundaries = self.build_inner_polygons()
    for innerBoundary in innerBoundaries:
        if outerBoundary.intersects_path(innerBoundary):
            innerBoundary = innerBoundary.clip_to_bbox((np.min(xt), np.min(yt), np.max(xt), np.max(yt)))
            mask = np.logical_or(zt.mask, innerBoundary.contains_points(xyt))
            zt = np.ma.masked_array(zt.data, mask)

    #TODO: Compound mask for other boundaries.
    if self.weir_boundaries is not None:
        for boundary in self.weir_boundaries:
            pass
    
    if self.culvert_boundaries is not None:
        for boundary in self.culvert_boundaries:
            pass

    if padding is not None:
        padding = griddata((padding[:,0], padding[:,1]), padding[:,2], (xt, yt), method='nearest')
        idx, = np.where(zt.mask)
        zt = np.ma.filled(zt, np.nan)
        zt[idx] = padding[idx]
    
    zt = zt.reshape(shape)
    return DEM.DEM(x, y, zt, geoTransform, self.epsg, self.datum)

def get_raster_from_extent(self, extent , dx, dy, epsg, padding=None):
    """
    This script rasterizes the mesh into a regular grid that can be exported as GeoTiff.
    Uselful for exploring the data on a GIS program.
    """
    min_x, max_x, min_y, max_y = extent
    
    x = np.arange(min_x, max_x+dx, dy)
    y = np.arange(min_y, max_y+dy, dy)
    
    xt, yt = np.meshgrid(x, y)
    shape = xt.shape
    
    xt = xt.reshape(xt.size)
    yt = np.flipud(yt.reshape(yt.size))
    xyt = np.vstack((xt,yt)).T

    bbox_path = Path([(np.min(xt), np.min(yt)),
                    (np.max(xt), np.min(yt)),
                    (np.max(xt), np.max(yt)),
                    (np.min(xt), np.max(yt)),
                    (np.min(xt), np.min(yt))], closed=True)

    xyz = self.get_xyz(extent=bbox_path, radius=0.2)
    zt = griddata((xyz[:,0], xyz[:,1]), xyz[:,2], (xt, yt), method='linear', fill_value=np.nan)
    
    
    # generate masks
    if ~hasattr(self,'outer_boundary'):
        self.outer_boundary = self.build_outer_polygon()
    path = self.outer_boundary.clip_to_bbox((np.min(xt), np.min(yt), np.max(xt), np.max(yt)))
    mask1 = np.isin(zt, -99999.0)
    mask2 = path.contains_points(xyt)
    mask = np.logical_or(mask1, ~mask2)
    for i in range(len(self.land_boundaries)):
        if self.land_boundaries[i][-1] in [1,21]:
            extents = self.land_boundaries[i][0].get_extents()
            extents = extents.get_points()
            c1 = np.logical_and(np.min(xt) <= extents[1,0], np.max(xt) >= extents[0,0])
            c2 = np.logical_and(np.min(yt) <= extents[1,1], np.max(yt) >= extents[0,1])
            if np.logical_and(c1,c2):
                shape = self.land_boundaries[i][0].clip_to_bbox((np.min(xt), np.min(yt), np.max(xt), np.max(yt)))
                mask3 = shape.contains_points(xyt)
                mask = np.logical_or(mask, mask3)  
    geoTransform = (np.min(min_x), dx, 0, np.max(y), 0, -np.abs(dy))
    depth = np.ma.masked_array(zt, mask).reshape(_shape)
    if padding is not None:
        depth = adcpy.adcirc.raster._apply_padding(x, y, depth, padding)
    DEM = DEM.DEM()
    return DEM(x, y, depth, geoTransform, 4326, self.datum)
