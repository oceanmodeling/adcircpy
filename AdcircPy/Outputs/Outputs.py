import os
import numpy as np
import fnmatch
from netCDF4 import Dataset
import numpy as np
from AdcircPy.Model import AdcircMesh
from AdcircPy.Outputs.Maxele import Maxele
from AdcircPy.Outputs.ElevationStations import ElevationStations
from AdcircPy.Outputs.ScalarSurfaceExtrema import ScalarSurfaceExtrema

class Outputs(object):
  """
  Private class called by AdcircPy.read_output() that returns
  the appropriate subclass belonging to the given output file.
  Supports ASCII and NetCDF outputs.
  Fortran binary outputs are not supported.
  """
  def __init__(self, path, fort14=None, datum='MSL', epsg=None, datum_grid=None):
    self.path = path
    self.datum=datum
    self.epsg=epsg
    self.datum_grid=datum_grid
    self.fort14=fort14
    self._init_fort14()

  def get_output_instance(self):
    if self._is_ncfile() == True:
      return self._netcdf_factory()
    else:
      return self._ascii_factory()

  def _init_fort14(self):
    if self.fort14 is None and self._is_ncfile()==False:
        raise Exception('For reading ASCII outputs, a fort.14 is required. A fort.15 is optional.')
    elif isinstance(self.fort14, str):
      self.fort14 = AdcircMesh.from_fort14(fort14=self.fort14, datum=self.datum, epsg=self.epsg, datum_grid=self.datum_grid)

  def _is_ncfile(self):
    try:
      Dataset(self.path)
      return True
    except:
      return False

  def _netcdf_factory(self):
    nc = Dataset(self.path)
    if 'station' in nc.variables.keys():
      if 'zeta' in nc.variables.keys():
        return ElevationStations.from_netcdf(self.path)
      else:
        raise NotImplementedError('The OutputFactory class has not implemented this output type yet, or this is not an Adcirc output file.')  
    elif 'zeta_max' in nc.variables.keys():
      return Maxele.from_netcdf(self.path, self.fort14, self.datum, self.epsg, self.datum_grid)
    elif 'zeta' in nc.variables.keys():
      return ElevationStations.from_netcdf(self.path)
    else:
      raise NotImplementedError('Guessed a NetCDF output but instantiation has not been implemented yet.')

  def _ascii_factory(self):
    self.f = open(self.path, 'r')
    # Might be a harmonic constituents file...
    if self.check_is_HarmonicConstituentsSurface()==True:
      return self._init_HarmonicConstituentsSurface()
    
    # else, could be a general surface file...
    elif self.check_is_SurfaceFile()==True:
      if self.number_of_points == self.fort14.x.size:
        if self.number_of_datasets > 2:
          return self._init_ScalarSurfaceTimeseries()
        elif self.number_of_datasets in [1,2]:
          return self._init_ScalarSurfaceMaxima()
    
    # else has to be a stations output file.
    else:
      return self._init_OutputStations()

  def check_is_HarmonicConstituentsSurface(self):
    self.line = self.f.readline().strip()
    try:
      int(self.line)
      return True
    except:
      return False

  def check_is_SurfaceFile(self):
    self.line = self.f.readline().split()
    self.number_of_datasets = int(self.line[0])
    self.number_of_points = int(self.line[1])
    if self.number_of_points==self.fort14.x.size:
      return True

  def _init_ScalarSurfaceTimeseries(self):
    """ """
    raise NotImplementedError("Guessed Scalar Surface Timeseries output, but instantiantion is not yet implemented")

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

  def _init_HarmonicConstituentsSurface(self):
    raise NotImplementedError('Guessed a Harmonic Constituents Output Surface but instantiantion has not yet been implemented.')

  def _init_OutputStations(self):
    raise NotImplementedError('Guessed an Output Stations file but instantiantion has not yet been implemented.')

  def __del__(self):
    if hasattr(self, 'f'):
      self.f.close()