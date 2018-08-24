import os
import numpy as np
import fnmatch
from netCDF4 import Dataset
import numpy as np
from AdcircPy.Mesh import AdcircMesh
from AdcircPy import Outputs


def read_outputs(path, **kwargs):
    return Outputs.Outputs(path, **kwargs)._open_file()

def _open_file(self):
  if os.path.isfile(self._path)==False:
    raise FileNotFoundError("No such file or directory: %s" % path)
  try:
    Dataset(self._path)
    _nc = True
  except:
    _nc = False
  if _nc == True:
    return self._read_netcdf()
  else:
    if self.fort14 is None:
      raise Exception('For reading ASCII outputs, a fort.14 is required. A fort.15 is optional.')
    if isinstance(self.fort14, AdcircMesh)==False:
      self.fort14 = AdcircMesh.from_fort14(fort14=self.fort14, fort15=self.fort15, datum=self.datum, epsg=self.epsg)
    return self._read_ascii()

def _read_netcdf(self):
  self.Dataset = Dataset(self._path)
  if 'station' in self.Dataset.dimensions.keys():
    if 'zeta' in self.Dataset.variables.keys():
      return Outputs.ElevationStations.from_netcdf(self._path)

def _read_ascii(self):
  f = open(self._path)
  line = f.readline().strip()
  try: int(line); _=True
  except: pass; _=False
  if _==True:
    f.close()  
    return self._harmonic_constituent_ascii()
  line = f.readline().split()
  number_of_datasets = int(line[0])
  number_of_points = int(line[1])
  if number_of_points == self.fort14.x.size:
    f.close()
    return Outputs.__OutputSurface.from_ascii(self._path, self.fort14, epsg=self.epsg)
  else:
    f.close()
    return Outputs.__OutputStations.from_ascii(self._path, self.fort14)

def _harmonic_constituent_ascii():
  pass
  



  # Error checking for input args
  # if fort14 is None:
  #     raise IOError("A fort.14 file is required to parse ASCII outputs.")
  # if isinstance(fort14, ("".__class__, u"".__class__)):
  #     fort14 = Mesh.init_from_fort14(fort14, datum, epsg)
  # elif isinstance(fort14, Mesh):
  #     pass
  # else:
  #     raise IOError("fort14 keyword provided is neither a path to a fort.14 ASCII file nor an AdcircPy.Mesh instance!")
  # _shape = fort14.values.shape
  # fort14 = fort14.get_dict()
  
  # f = open(path)
  # description          = f.readline().strip()
  # line                 = f.readline().split()
  # number_of_datasets   = int(line[0]) 
  # number_of_datapoints = int(line[1])# This is either NP or Number of stations
  # output_time          = float(line[2])
  # output_interval      = float(line[3])
  # record_type          = int(line[4]) # 1 for elevation, 2 for velocity, 3 for 3D
  
  # # This is probably a gridded output
  # if number_of_datapoints==_shape[0] and number_of_datasets>2:
  #     time     = list()
  #     timestep = list()
  #     values   = list()
  #     nodeID   = list()
  #     for i in range(number_of_datasets):
  #         line = f.readline().split()
  #         time.append(float(line[0]))
  #         timestep.append(int(line[1]))
  #         for i in range(number_of_datapoints):
  #             _values = list()
  #             line = f.readline().split()
  #             nodeID.append(int(line[0]))
  #             if record_type == 1:
  #                 _values.append(float(line[1]))
  #             elif record_type == 2:
  #                 _values.append((float(line[1]),float(line[2])))
  #             elif record_type == 3:
  #                 line = f.readline().split()
  #                 _values.append((float(line[1]), float(line[2]), float(line[3])))
  #         _values = np.asarray(_values)
  #         _values = np.ma.masked_equal(_values, -99999.0)
  #         values.append(_values)
  #     f.close()
  #     fort14['values'] = values
  #     return Outputs.SurfaceTimeseries(**fort14)
  
  # # this is probably a station timeseries
  # elif number_of_datapoints<_shape[0] and number_of_datasets>2:
  #     stations = dict()
  #     time=list()
  #     timestep=list()
  #     for i in range(number_of_datasets):
  #         line = f.readline().split()
  #         time.append(float(line[0]))
  #         timestep.append(int(line[1]))
  #         for j in range(number_of_datapoints):
  #             if j not in stations.keys():
  #                 stations[j]=list()
  #             stations[j].append(float(f.readline().split()[-1]))
  #     f.close()
  #     return Outputs.StationTimeseries()
  
  # # this is probably an extrema file (*.63)
  # elif number_of_datasets==1: 
  #     nodeID = list()
  #     values = list()
  #     for i in range(number_of_datasets):
  #         line = f.readline().split()
  #         for i in range(number_of_datapoints):
  #             line = f.readline().split()
  #             nodeID.append(int(line[0]))
  #             if record_type == 1:
  #                 values.append(float(line[1]))
  #             elif record_type == 2:
  #                 values.append((float(line[1]),float(line[2])))
  #             elif record_type == 3:
  #                 values.append((float(line[1]), float(line[2]), float(line[3])))
  #         values = np.asarray(values)
  #         values = np.ma.masked_equal(values, -99999.0)
  #     fort14['values'] = values
  #     f.close()
  #     return Outputs.SurfaceExtrema(**fort14)
  
  # elif number_of_datasets==2:
  #     pass

        

  #   