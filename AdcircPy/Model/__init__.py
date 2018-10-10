from collections import OrderedDict
from scipy.spatial import KDTree
from AdcircPy.Model import UnstructuredMesh as _UnstructuredMesh
from AdcircPy.Model import AdcircMesh as _AdcircMesh
from AdcircPy.Model import AdcircRun as _AdcircRun
from AdcircPy.Model import ElevationStationsOutput as _ElevationStationsOutput
from AdcircPy.Model import _fort14
from AdcircPy.Model import _fort13

class UnstructuredMesh(object):   
  def __init__(self, x, y, elements, values,
                                     nodeID=None,
                                     elementID=None,
                                     datum=None,
                                     epsg=None,
                                     datum_grid=None,
                                     ocean_boundaries=None,
                                     land_boundaries=None,
                                     inner_boundaries=None,
                                     weir_boundaries=None,
                                     inflow_boundaries=None,
                                     outflow_boundaries=None,
                                     culvert_boundaries=None):
    self._x           = x
    self._y           = y
    self._values      = values
    self.elements     = elements
    self.nodeID       = nodeID
    self.elementID    = elementID
    self.datum        = datum
    self.epsg         = epsg
    self._datum_grid  = datum_grid
    self.ocean_boundaries    = ocean_boundaries
    self.land_boundaries     = land_boundaries
    self.inner_boundaries    = inner_boundaries
    self.weir_boundaries     = weir_boundaries
    self.inflow_boundaries   = inflow_boundaries
    self.outflow_boundaries  = outflow_boundaries
    self.culvert_boundaries  = culvert_boundaries
    self._init_datum_grid()
    self._init_Tri()
    self._init_KDTree()

  @property
  def values(self):
    return self._values

  @property
  def x(self):
    return self._x
  
  @property
  def y(self):
    return self._y
    
  def get_mean_value(self, **kwargs):
    """   """
    return _UnstructuredMesh.get_mean_value(self, **kwargs)

  def get_contours(self, levels, **kwargs):
    """   """
    return _UnstructuredMesh.get_contours(self, levels, **kwargs)

  def get_values_under_Path(self, Path, **kwargs):
    """
    Input is a matplotlib.path.Path instance.
    Returns array of values under specified path.
    """
    return _UnstructuredMesh.get_values_under_Path(self, Path, **kwargs)

  def make_plot(self, **kwargs):
    """  """
    return _UnstructuredMesh.make_plot(self, **kwargs)
  
  def get_values_at_xy(self, x, y, **kwargs):
    """
    Returns numpy array of values at coordinates or list of coordinates give by lon lat
    """
    return _UnstructuredMesh.get_values_at_xy(self, x, y, **kwargs)

  def rasterize_to_geoTransform(self, geoTransform, shape, **kwargs):
    """  """
    return _UnstructuredMesh.rasterize_to_geoTransform(self, geoTransform, shape, **kwargs)

  def get_xyz(self, extent=None, **kwargs):
    """
    Returns numpy array of dimension [D,3] with columns <x, y, values>
    """
    return _UnstructuredMesh.get_xyz(self, **kwargs)
  
  def get_xy(self, extent=None, epsg=None):
    """
    Returns numpy array of dimension [D,2] with columns <x, y>
    """
    return _UnstructuredMesh.get_xy(self, extent, epsg)

  def get_elements_from_exent(self, extent):
    """ """
    return _UnstructuredMesh.get_elements_from_extent(self, extent)

  def get_finite_volume(self, index):
    """ """
    return _UnstructuredMesh.get_finite_volume(self, index)
    
  def get_Path_from_element(self, element):
    """ """
    return _UnstructuredMesh.get_Path_from_element(self, element)
  
  def get_extent(self, **kwargs):
    """ """
    return _UnstructuredMesh.get_extent(self, **kwargs)

  def get_element_containing_coord(self, coord):
    return _UnstructuredMesh.get_element_containing_coord(self, coord)

  def get_finite_volume_Path_list(self, index):
    return _UnstructuredMesh.get_finite_volume_Path_list(self, index)

  def get_finite_volume_element_list(self, index):
    return _UnstructuredMesh.get_finite_volume_element_list(self, index)
  
  def transform_to_epsg(self, epsg):
    """ """
    _UnstructuredMesh.transform_to_epsg(self, epsg)

  def plot_trimesh(self, **kwargs):
    """ """
    return _UnstructuredMesh.plot_trimesh(self, **kwargs) 
  
  def plot_outerBoundary(self):
    return _UnstructuredMesh.plot_outerBoundary(self)
  
  def _init_fig(self, axes=None, extent=None, title=None, epsg=None):
    """    """
    _UnstructuredMesh._init_fig(self, axes, extent, title, epsg)

  def _init_cbar(self, cmap, vmin, vmax):
    """ """
    _UnstructuredMesh._init_cbar(self, cmap, vmin, vmax)

  def _init_Tri(self):
    _UnstructuredMesh._init_Tri(self)

  def _init_KDTree(self):
    self.KDTree = KDTree(self.get_xy())

  def _init_datum_grid(self):
    _UnstructuredMesh._init_datum_grid(self)

  def _get_finite_volume_interp(self, point_idx):
    """
    Fetches the Path for the DEM interpolator.
    """
    return _UnstructuredMesh._get_finite_volume_interp(self, point_idx)

  def _get_extent_idx(self, extent, epsg, **kwargs):
    """
    Finds the indices of the mesh nodes inside a bounding box.
    kwargs:
        extent : list of the form [min_x, max_x, min_y, max_y]
    """
    return _UnstructuredMesh._get_extent_idx(self, extent, epsg, **kwargs)

  def __sub__(self, other):
    """
    Used in the case where two different grids are subtracted directly.
    Grids must have the same dimensions.
    Example:
        grid1 = adcpy.read_grid("fort.14.1")
        grid2 = adcpy.read_grid("fort.14.2")
        diff = grid2 - grid1
        diff.make_plot(show=True)
    """
    return _UnstructuredMesh.get_difference(self, other)

class AdcircMesh(UnstructuredMesh):
  def __init__(self, x, y, elements, values,
                                     fort13=None,
                                     fort15=None,
                                     nodeID=None,
                                     description=None,
                                     elementID=None,
                                     datum=None,
                                     epsg=None,
                                     datum_grid=None,
                                     ocean_boundaries=None,
                                     land_boundaries=None,
                                     inner_boundaries=None,
                                     weir_boundaries=None,
                                     inflow_boundaries=None,
                                     outflow_boundaries=None,
                                     culvert_boundaries=None):
    UnstructuredMesh.__init__(self, x, y, elements, values, nodeID, elementID,
                              datum, epsg, datum_grid, ocean_boundaries,
                              land_boundaries, inner_boundaries, weir_boundaries,
                              inflow_boundaries, outflow_boundaries, culvert_boundaries)
    self.description = description
    self._init_fort13(fort13)

  @classmethod
  def from_fort14(cls, fort14, datum='MSL', epsg=4326, fort13=None, fort15=None, datum_grid=None):
    return _AdcircMesh.from_fort14(cls, fort14, datum, epsg, fort13, fort15, datum_grid)

  def make_plot(self, **kwargs):
    return _AdcircMesh.make_plot(self, **kwargs)

  def interpolate_DEM(self, DEM, **kwargs):
    _AdcircMesh.interpolate_DEM(self, DEM, **kwargs)
  
  def write_fort14(self, path):
    _fort14.write_fort14(self, path)

  def _init_fort13(self, path):
    _fort13._init_fort13(self, path)

class AdcircRun(object):

  def __init__(self, AdcircMesh, Tides=None, Winds=None, Waves=None,
                ElevationStationsOutput=None, VelocityStationsOutput=None,
                ElevationGlobalOutput=None, VelocityGlobalOutput=None, 
                **kwargs):
    self.AdcircMesh = AdcircMesh
    self.init_Tides(Tides)
    self.init_Winds(Winds)
    # self.init_Waves(Waves)
    self.ElevationStationsOutput = ElevationStationsOutput
    self.init_fort15(**kwargs)

  def init_Tides(self, Tides):
    _AdcircRun.init_Tides(self, Tides)

  def init_Winds(self, Winds):
    self.Winds = Winds

  def init_fort15(self):
    _AdcircRun.init_fort15(self)

  def dump(self, directory):
    _AdcircRun.dump(self, directory)

  def _init_TPXO(self):
    _AdcircRun._init_TPXO(self)

  def _write_fort15(self):
    _AdcircRun._write_fort15(self)

  def _write_IHOT(self):
    _AdcircRun._write_IHOT(self)

  def _write_NWP(self):
    _AdcircRun._write_NWP(self)

  def _write_NWS(self):
    _AdcircRun._write_NWS(self)

  def _write_NRAMP(self):
    _AdcircRun._write_NRAMP(self)

  def _write_TAU0(self):
    _AdcircRun._write_TAU0(self)

  def _write_DTDP(self):
    _AdcircRun._write_DTDP(self)

  def _write_RNDAY(self):
    _AdcircRun._write_RNDAY(self)

  def _write_DRAMP(self):
    _AdcircRun._write_DRAMP(self)

  def _write_H0_VELMIN(self):
    _AdcircRun._write_H0_VELMIN(self)

  def _write_SLAM0_SFEA0(self):
    _AdcircRun._write_SLAM0_SFEA0(self)

  def _write_FFACTOR(self):
    _AdcircRun._write_FFACTOR(self)

  def _write_ESLM(self):
    _AdcircRun._write_ESLM(self)

  def _write_CORI(self):
    _AdcircRun._write_CORI(self)

  def _write_NTIF(self):
    _AdcircRun._write_NTIF(self)

  def _write_NBFR(self):
    _AdcircRun._write_NBFR(self)

  def _write_station_outputs(self):
    _AdcircRun._write_station_outputs(self)

  def _write_global_outputs(self):
    _AdcircRun._write_global_outputs(self)

  def _write_harmonic_outputs(self):
    _AdcircRun._write_harmonic_outputs(self)

  def _write_hotstart_parameters(self):
    _AdcircRun._write_hotstart_parameters(self)

  def _write_iteration_parameters(self):
    _AdcircRun._write_iteration_parameters(self)

  def _write_netcdf_parameters(self):
    _AdcircRun._write_netcdf_parameters(self)

  def _write_fortran_namelists(self):
    _AdcircRun._write_fortran_namelists(self)

class _AdcircOutputs(object):
 def __init__(self, stations, netcdf):
    self.stations = stations
    self.netcdf = netcdf

class ElevationStationsOutput(_AdcircOutputs):
  def __init__(self, stations, netcdf=True):
    super(ElevationStationsOutput, self).__init__(stations, netcdf)
    self.init_params()

  @classmethod
  def from_fort15(self, path):
    _ElevationStationsOutput.from_fort15(self, path)

  @classmethod
  def from_csv(self, path):
    _ElevationStationOutput.from_csv(self, path)  

  def init_params(self):
    _ElevationStationOutput.init_params(self)






# class VelocityStationOutputs(dict):
  
#   def __init__(self, netcdf=True):
#     self.netcdf=netcdf
  
#   @classmethod
#   def from_fort15(self, path):
#     _ElevationStationOutputs.from_fort15(self, path)

#   @classmethod
#   def from_csv(self, path):
#     _ElevationStationOutputs.from_csv(self, path)

# class ElevationGlobalOutputs(dict):

#   def __init__(self, netcdf=True):
#     self.netcdf=netcdf
  
#   @classmethod
#   def from_fort15(self, path):
#     _ElevationStationOutputs.from_fort15(self, path)

#   @classmethod
#   def from_csv(self, path):
#     _ElevationStationOutputs.from_csv(self, path)

# class VelocityGlobalOutputs(dict):
  
#   def __init__(self, netcdf=True):
#     self.netcdf=netcdf
  
#   @classmethod
#   def from_fort15(self, path):
#     _ElevationStationOutputs.from_fort15(self, path)

#   @classmethod
#   def from_csv(self, path):
#     _ElevationStationOutputs.from_csv(self, path)