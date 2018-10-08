import os
import numpy as np
import fnmatch
from netCDF4 import Dataset
import numpy as np
from AdcircPy.Model import AdcircMesh
from AdcircPy import Outputs

def __new__(cls, path, fort14=None, fort15=None, datum='MSL', epsg=None, datum_grid=None):

  # Read and return NetCDF output
  if cls.is_ncfile(path) == True:
    return cls.netcdf_factory(path, fort14, fort15, datum, epsg, datum_grid)

  # Read and return ASCII output
  else:
    if fort14 is None:
      raise Exception('For reading ASCII outputs, a fort.14 is required. A fort.15 is optional.')
    if isinstance(fort14, AdcircMesh)==False:
      fort14 = AdcircMesh.from_fort14(fort14=fort14, fort15=fort15, datum=datum, epsg=epsg, datum_grid=datum_grid)
    f = open(path)
    line = f.readline().strip()
    
    # Test to see if it's a harmonic constituents file
    try: int(line); _=True
    except: pass; _=False
    if _==True:
      f.close()  
      return Outputs.HarmonicConstituentsSurface()
    
    # Test which other types this file might be
    line = f.readline().split()
    number_of_datasets = int(line[0])
    number_of_points = int(line[1])
    # Correspods to a surface file
    if number_of_points == fort14.x.size:
      # Correspondos to a surface time series file
      if number_of_datasets > 2:
        return Outputs._ScalarSurfaceTimeseries(fort14.x,
                                                     fort14.y,
                                                     values,
                                                     fort14.elements,
                                                     f,
                                                     epsg=epsg)

      # Corresponds to a surface maxima file
      elif number_of_datasets in [1,2]:
        line = f.readline().split()
        _time = float(line[0].strip(' /n'))
        _timestep = int(line[1])
        _nodeid = list()
        _values = list()
        for i in range(number_of_points):
          line = f.readline().split()
          _nodeid.append(int(line[0].strip(' \n')))
          _values.append(float(line[1].strip(' \n')))
        _extrema_time=list()
        if number_of_datasets==2:
          for i in range(number_of_points):
            line = f.readline().split()
            _extrema_time.append(float(line[1].strip(' \n')))
        nodeID = np.asarray(_nodeid)
        values = np.ma.masked_equal(_values, -99999.)
        f.close()
        return Outputs._ScalarSurfaceExtrema(fort14.x,
                                             fort14.y,
                                             fort14.elements,
                                             values,
                                             _extrema_time,
                                             epsg=epsg,
                                             nodeID=nodeID)
    # Then it has to be a station timeseries file
    else:
      f.close()
      return Outputs._OutputStations.from_ascii(self._path, self._fort14)

def is_ncfile(path):
  try: Dataset(path); nc = True
  except: nc = False
  return nc

def netcdf_factory(path, fort14=None, fort15=None, datum='MSL', epsg=None, datum_grid=None):
  nc = Dataset(path)
  if 'station' in nc.variables.keys():
    if 'zeta' in nc.variables.keys():
      return ElevationStations.from_netcdf(path)
    else:
      raise NotImplementedError('The OutputFactory class has not implemented this output type yet, or this is not an Adcirc output file.')  
  elif 'zeta_max' in nc.variables.keys():
    return Outputs.Maxele.from_netcdf(path, fort14, fort15, datum, epsg, datum_grid)
  else:
    raise NotImplementedError('The OutputFactory class has not implemented this output type yet, or this is not an Adcirc output file.')

def ascii_factory(path, fort14=None, fort15=None, datum='MSL', epsg=None, datum_grid=None):
  pass