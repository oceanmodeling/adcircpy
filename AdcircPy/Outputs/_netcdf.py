from netCDF4 import Dataset
from AdcircPy import Outputs

def _set__nc(self):
    try:
        Dataset(self._path)
        self._nc = True
    except:
        self._nc = False

def _get_output(self):
    if 'station' in self.Dataset.dimensions.keys():
        if 'zeta' in self.Dataset.variables.keys():
            return Outputs.ElevationStationTimeSeries.from_netcdf(self._path)

def _read_netcdf(self):
    self.Dataset = Dataset(self._path)
    return self._get_netcdf_output()

def _add_fort14_data(self, params):
    pass
























    # if 'station_name' in nc.variables.keys():
    #     stations = dict()     
    #     for i, station_id in enumerate(nc.variables['station_name'][:]):
    #         id=''.join(station_id).strip()
    #         stations[id]=dict()
    #         stations[id]['coordinates'] = (nc.variables['x'][i], nc.variables['y'][i])
    #         stations[id]['values'] = nc.variables['zeta'][:,i]
    #         stations[id]['time'] = nc.variables['time'][:]
    #     return adcpy.adcirc.timeseries(**stations)

    # elif 'adcirc_mesh' in nc.variables.keys():

    #     if fort14 is not None:
    #         if isinstance(fort14, ("".__class__, u"".__class__)):
    #             fort14 = adcpy.read_mesh(fort14, fort13, fort15, datum, epsg)
            
    #         elif type(fort14)==type(adcpy.adcirc.topobathy()):
    #             params = {  'x'                  : fort14.x,
    #                         'y'                  : fort14.y,
    #                         'elements'           : fort14.elements,
    #                         'values'             : fort14.values,
    #                         'nodeID'             : fort14.nodeID,
    #                         'elementID'          : fort14.elementID, 
    #                         "oceanBoundaries"    : fort14.oceanBoundaries,
    #                         "landBoundaries"     : fort14.landBoundaries,
    #                         "innerBoundaries"    : fort14.innerBoundaries,
    #                         "weirBoundaries"     : fort14.weirBoundaries,
    #                         "inflowBoundaries"   : fort14.inflowBoundaries,
    #                         "outflowBoundaries"  : fort14.outflowBoundaries,
    #                         "culvertBoundaries"  : fort14.culvertBoundaries,
    #                         "fort14_description" : fort14.fort14_description,
    #                         "fort14_path"        : fort14.fort14_path}
    #         else:
    #             raise IOError("fort14 keyword provided is neither a path to a fort.14 file nor an adcirc topobathy instance!")

    #     else:
    #         # Rebuild mesh from ncoutput if mesh file is not provided.
    #         ocean_boundaries = list()
    #         _ocean_boundaries = nc.variables['nbdv'][:]
    #         for i in range(_ocean_boundaries.shape[1]):
    #             ocean_boundaries.append(nc.variables['nbdv'][:][:,i].flatten())

    #         land_boundaries    = list()
    #         inner_boundaries   = list()
    #         inflow_boundaries  = list()
    #         outflow_boundaries = list()
    #         weir_boundaries    = list()
    #         culvert_boundaries = list()

    #         for i, btype in enumerate(nc.variables['ibtype'][:]):
    #             node_index = nc.variables['nbvv'][:][:,i] - 1
    #             node_index = list(filter(lambda x: x!=-1, node_index))
    #             if btype in [0, 10, 20]:
    #                 land_boundaries.append([np.asarray(node_index).flatten(), btype])
    #             elif btype in [1, 11, 21]:
    #                 inner_boundaries.append([np.asarray(node_index).flatten(), btype])
    #             elif btype in [2, 12, 22, 102, 122]:
    #                 inflow_boundaries.append([np.asarray(node_index).flatten(), btype])
    #             elif btype in [3, 13, 23]:
    #                 outflow_boundaries.append([node_index, btype])
    #             elif btype in [4, 24]:
    #                 _weir_boundaries = dict()
    #                 _weir_boundaries['front_face'] = np.asarray(node_index).flatten()
    #                 _weir_boundaries['back_face']  = None
    #                 _weir_boundaries['height']     = None
    #                 _weir_boundaries['subcritical_flow_coefficient'] = None
    #                 _weir_boundaries['supercritical_flow_coefficient'] = None
    #                 _weir_boundaries['btype'] = btype
    #                 weir_boundaries.append(_weir_boundaries)
    #             elif btype in [5, 25]:
    #                 _culvert_boundaries = dict()
    #                 _culvert_boundaries['front_face'] = np.asarray(node_index).flatten()
    #                 _culvert_boundaries['back_face']  = None
    #                 _culvert_boundaries['height'] = None
    #                 _culvert_boundaries['subcritical_flow_coefficient'] = None
    #                 _culvert_boundaries['supercritical_flow_coefficient'] = None
    #                 _culvert_boundaries['cross_barrier_pipe_height'] = None
    #                 _culvert_boundaries['friction_factor'] = None
    #                 _culvert_boundaries['pipe_diameter'] = None
    #                 _culvert_boundaries['btype'] = btype
    #                 culvert_boundaries.append(_culvert_boundaries)
    
    #         params = {  'x'                 : nc.variables['x'][:].flatten(),
    #                     'y'                 : nc.variables['y'][:].flatten(),
    #                     'elements'          : nc.variables['element'][:]-1,
    #                     'nodeID'            : np.arange(nc.variables['x'][:].shape[0]).flatten(),
    #                     'elementID'         : np.arange(nc.variables['element'][:][:,0].shape[0]).flatten(),
    #                     'time'              : nc.variables['time'][:].flatten(),
    #                     'oceanBoundaries'   : ocean_boundaries,
    #                     'landBoundaries'    : land_boundaries,
    #                     'innerBoundaries'   : inner_boundaries,
    #                     'inflowBoundaries'  : inflow_boundaries,
    #                     'outflowBoundaries' : outflow_boundaries,
    #                     'weirBoundaries'    : weir_boundaries,
    #                     'culvertBoundaries' : culvert_boundaries}

    #     params['values'] = list()
    #     value_key=list()
    #     for key in nc.variables.keys():
    #         if key not in ['x', 'y', 'element', 'adcirc_mesh', 'neta',
    #                        'nvdll', 'max_nvdll', 'ibtypee','nbdv', 'nvel',
    #                        'nvell', 'max_nvell','ibtype','nbvv', 'depth',  'time']:
    #             value_key.append(key)
        
    #     if 'zeta_max':
    #         time_of_peak_key = fnmatch.filter(value_key, 'time_of_*').pop()
    #         params['values'] = nc.variables['zeta_max'][:]
    #         params['units']  = nc.variables['zeta_max'].units
    #         return Outputs.Maxele(**params)

    #     elif 'u_vel' in value_key:
    #         params['u'] = nc.variables['u_vel'][:]
    #         params['v'] = nc.variables['v_vel'][:]
    #         return Outputs.Velocity(**params)

    #     else:
    #         raise NotImplementedError("Definition for this type of NetCDF not found. Please add definition to _netcdf.py")
