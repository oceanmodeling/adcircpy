import os
import numpy as np
import fnmatch
from netCDF4 import Dataset
import numpy as np
from AdcircPy.Model import AdcircMesh
from AdcircPy.Outputs.Maxele import Maxele
from AdcircPy.Outputs.ElevationSurfaceTimeseries import ElevationSurfaceTimeseries
from AdcircPy.Outputs.ElevationStations import ElevationStations
from AdcircPy.Outputs.ScalarSurfaceExtrema import ScalarSurfaceExtrema
from AdcircPy.Outputs.HarmonicConstituentsElevationStations import HarmonicConstituentsElevationStations

class _OutputFactory(object):
  """
  Private class called by AdcircPy.read_output() that returns
  the appropriate subclass belonging to the given output file.
  Supports ASCII and NetCDF outputs.
  Fortran binary outputs are not supported.
  """
  def __init__(self, path, fort14=None, datum='MSL', epsg=None, datum_grid=None, fort15=None):
    """
    Should probably be replaced with __call__ in order to avoid having to call get_output_instance()
    """
    self.path = path
    self.datum=datum
    self.epsg=epsg
    self.datum_grid=datum_grid
    self.fort14=fort14
    self.fort15=fort15
    self._init_fort14()

  def get_output_instance(self):
    if self._is_ncfile() == True:
      return self._netcdf_factory()
    else:
      return self._ascii_factory()

  def _init_fort14(self):
    if isinstance(self.fort14, str):
      self.fort14 = AdcircMesh.from_fort14(fort14=self.fort14, datum=self.datum, epsg=self.epsg, datum_grid=self.datum_grid)

  def _is_ncfile(self):
    try:
      Dataset(self.path)
      return True
    except:
      return False

  def _netcdf_factory(self):
    nc = Dataset(self.path)
    if 'adcirc_mesh' in nc.variables.keys():
      if 'zeta_max' in nc.variables.keys():
        return Maxele.from_netcdf(self.path, self.fort14, self.datum, self.epsg, self.datum_grid)
      elif 'zeta' in nc.variables.keys():
        return ElevationSurfaceTimeseries.from_netcdf(self.path, self.fort14, self.datum, self.epsg, self.datum_grid)
      else:
        raise NotImplementedError('The OutputFactory class has not implemented this output type yet, or this is not an Adcirc output file.')  

    else:
      if 'zeta' in nc.variables.keys():
        return ElevationStations.from_netcdf(self.path)
      elif 'phs' in nc.variables.keys():
        return HarmonicConstituentsElevationStations.from_netcdf(self.path)
      else:
        raise NotImplementedError('The OutputFactory class has not implemented this output type yet, or this is not an Adcirc output file.')  

  def _ascii_factory(self):
    self.f = open(self.path, 'r')
    # Might be a harmonic constituents file...
    if self._check_is_HarmonicConstituentsFile()==True:
      number_of_points = self._check_is_SurfaceOrStations()
      if self.fort14 is not None and self.fort14.x.size==number_of_points:
        raise NotImplementedError('Guessed an ASCII HarmonicConstituentsOutputSurface but instantiantion has not yet been implemented.')
      elif self.fort15 is not None:
        # should probably parse fort.15 and compare the number of stations first.
        return HarmonicConstituentsElevationStations.from_ascii(self.path, self.fort15)
      else:
        raise Exception('When loading an ASCII Harmonic constituents file, the fort.14 needs to be provided for fort.53 files or the fort.15 needs to be provided for fort.51 files.')    
    
    # else, could be a general surface file...
    elif self._check_is_SurfaceFile()==True:
      if self.number_of_points == self.fort14.x.size:
        # Could be a surface time series
        if self.number_of_datasets > 2:
          raise NotImplementedError("Guessed Scalar Surface Timeseries output, but instantiantion is not yet implemented")
        # or a surface maxima file
        elif self.number_of_datasets in [1,2]:
          return self._init_ScalarSurfaceMaxima()
    
    # else has to be a stations output file.
    else:
      raise NotImplementedError('Guessed an Output Stations file but instantiantion has not yet been implemented.')

  def _check_is_HarmonicConstituentsFile(self):
    self.line = self.f.readline().strip()
    try:
      int(self.line)
      return True
    except:
      return False

  def _check_is_SurfaceOrStations(self):
    for i in range(int(self.line)):
      self.f.readline()
    return int(self.f.readline())

  def _check_is_SurfaceFile(self):
    # we keep pushing this check further and further below...
    if self.fort14 is None:
      raise Exception('A fort.14 is required for reading some output file types.')
    self.line = self.f.readline().split()
    self.number_of_datasets = int(self.line[0])
    self.number_of_points = int(self.line[1])
    if self.number_of_points==self.fort14.x.size:
      return True

  def _init_ScalarSurfaceMaxima(self):
    """ """
    self.line = self.f.readline().split()
    time = float(self.line[0].strip(' /n'))
    timestep = int(self.line[1])
    nodeID = list()
    values = list()
    for i in range(self.number_of_points):
      self.line = self.f.readline().split()
      nodeID.append(int(self.line[0].strip(' \n')))
      values.append(float(self.line[1].strip(' \n')))
    extrema_time_vector=list()
    if self.number_of_datasets==2:
      for i in range(self.number_of_points):
        self.line = self.f.readline().split()
        extrema_time_vector.append(float(self.line[1].strip(' \n')))
    nodeID = np.asarray(nodeID)
    values = np.ma.masked_equal(values, -99999.)
    return ScalarSurfaceExtrema(self.fort14.x,
                                 self.fort14.y,
                                 self.fort14.elements,
                                 values,
                                 extrema_time_vector,
                                 epsg=self.epsg,
                                 nodeID=nodeID)

  def __del__(self):
    if hasattr(self, 'f'):
      self.f.close()