import numpy as np
from AdcircPy.Mesh import AdcircMesh
from AdcircPy import Outputs

def _load_fort14(self):
  if isinstance(self.fort14, ("".__class__, u"".__class__)):
      self.fort14 = Mesh.init_from_fort14(fort14, datum, epsg)
  
  elif isinstance(fort14, Mesh):
      pass

def _read_ascii_type(self):
  f = open(self._path)
  line = f.readline().strip()
  try:
    _num = int(line)
    self._ascii_type = 'harmonic_constituents' 
  except:
    self._ascii_type = None

  if self._ascii_type == None:
    # detect other than harmonic constituents
    raise NotImplementedError('need to detect other than harmonic constituents')
  f.close()





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