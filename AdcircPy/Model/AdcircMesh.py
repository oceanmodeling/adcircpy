import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import griddata
from AdcircPy.core import FixPointNormalize
from AdcircPy.Model.UnstructuredMesh import UnstructuredMesh
from AdcircPy.Model.NodalAttributes import NodalAttributes
from AdcircPy.Model.TidalRun import TidalRun
from AdcircPy.Model.HindcastRun import HindcastRun

class AdcircMesh(UnstructuredMesh):
  def __init__(self, x, y, elements, values, fort13=None, description=None, **kwargs):
    self.description = description
    self._init_fort13(fort13, kwargs.pop("spinup_attributes", None), kwargs.pop("runtime_attributes", None))
    super(AdcircMesh, self).__init__(x, y, elements, values, **kwargs)

  @classmethod
  def from_fort14(cls, fort14, datum='MSL', epsg=4326, fort13=None, datum_grid=None):
    kwargs = cls.parse_fort14(fort14)
    kwargs['datum'] = datum
    kwargs['epsg']  = epsg
    kwargs['fort13'] = fort13
    kwargs['datum_grid'] = datum_grid
    return cls(**kwargs)

  def TidalRun(self, start_date, end_date, spinup_date=None, constituents=None, netcdf=True, **kwargs):
    """ Instantiates an ADCIRC tidal only run. """
    return TidalRun(self, start_date, end_date, spinup_date, constituents, netcdf, **kwargs)
  
  def HindcastRun(self, hurdat_id, start_date=None, end_date=None, spinup_date=None, tides=True, netcdf=True, **kwargs):
    """ Generates an ADCIRC hindcast run using wind data from the HURDAT2 database. """
    return HindcastRun(self, hurdat_id, start_date, end_date, spinup_date, tides, netcdf, **kwargs)

  @staticmethod
  def parse_fort14(path):
    f  = open(path) 
    description = f.readline().rstrip()
    NE, NP = map(int, f.readline().split())  
    nodeID = list()
    x      = list()
    y      = list()
    z      = list()
    for k in range(NP):
        line = f.readline().split()
        nodeID.append(int(line[0])-1)
        x.append(float(line[1]))
        y.append(float(line[2]))
        z.append(float(line[3]))
    nodeID = np.squeeze(nodeID)
    x    = np.squeeze(x)
    y    = np.squeeze(y)
    z  = (-1)*np.squeeze(z)
    elementID = list()
    elements  = list()
    for k in range(NE):
        line = f.readline().split()
        elementID.append(int(line[0]))
        elements.append([int(line[2])-1, int(line[3])-1, int(line[4])-1])
    elementID = np.squeeze(elementID)
    elements  = np.array(elements)
    try:
        number_of_ocean_boundaries = int(f.readline().split()[0])
        total_number_of_ocean_boundary_nodes = int(f.readline().split()[0])
        ocean_boundaries = list()
        for i in range(number_of_ocean_boundaries):
            Number_of_nodes_for_ocean_boundary_x = int(f.readline().split()[0])
            node_index=list()
            for j in range(Number_of_nodes_for_ocean_boundary_x):
                node_index.append(int(f.readline())-1)
            node_index = np.asarray(node_index).flatten()
            ocean_boundaries.append(node_index)

        number_of_land_boundaries = int(f.readline().split()[0])
        total_number_of_land_boundary_nodes = int(f.readline().split()[0])
        
        land_boundaries = list()
        inner_boundaries  = list()
        inflow_boundaries = list()
        outflow_boundaries = list()
        weir_boundaries = list()
        culvert_boundaries = list()
        for i in range(number_of_land_boundaries):
            line = f.readline().split()
            number_of_nodes_for_land_boundary_x = int(line[0])
            btype = int(line[1])
            node_index=list()
            weir0 = list()
            weir1 = list()
            weir2 = list()
            weir3 = list()
            weir4 = list()
            weir5 = list()
            weir6 = list()
            weir7 = list()
            for j in range(number_of_nodes_for_land_boundary_x):
                line = f.readline().split()
                if btype in [0, 1, 2, 10, 11, 12, 20, 21, 22, 30, 102, 122]:
                    node_index.append(int(line[0])-1)
                elif btype in [3, 13, 23]:
                    node_index.append([int(line[0])-1, -float(line[1]), float(line[2])])
                elif btype in [4, 24]:
                    weir0.append(int(line[0])-1)
                    weir1.append(int(line[1])-1)
                    weir2.append(float(line[2]))
                    weir3.append(float(line[3]))
                    weir4.append(float(line[4]))   
                elif btype in [5, 25]:
                    weir0.append(int(line[0])-1)
                    weir1.append(int(line[1])-1)
                    weir2.append(float(line[2]))
                    weir3.append(float(line[3]))
                    weir4.append(float(line[4]))
                    weir5.append(float(line[5]))
                    weir6.append(float(line[6]))
                    weir7.append(float(line[7]))
                    
            if btype in [0, 10, 20]:
                land_boundaries.append([np.asarray(node_index).flatten(), btype])
            elif btype in [1, 11, 21]:
                inner_boundaries.append([np.asarray(node_index).flatten(), btype])
            elif btype in [2, 12, 22, 102, 122]:
                inflow_boundaries.append([np.asarray(node_index).flatten(), btype])
            elif btype in [3, 13, 23]:
                outflow_boundaries.append([node_index, btype])
            elif btype in [4, 24]:
                _weir_boundaries = dict()
                _weir_boundaries['front_face'] = np.asarray(weir0).flatten()
                _weir_boundaries['back_face']  = np.asarray(weir1).flatten()
                _weir_boundaries['height']     = np.asarray(weir2).flatten()
                _weir_boundaries['subcritical_flow_coefficient'] = np.asarray(weir3).flatten()
                _weir_boundaries['supercritical_flow_coefficient'] = np.asarray(weir4).flatten()
                _weir_boundaries['btype'] = btype
                weir_boundaries.append(_weir_boundaries)
            elif btype in [5, 25]:
                _culvert_boundaries = dict()
                _culvert_boundaries['front_face'] = np.asarray(weir0).flatten()
                _culvert_boundaries['back_face']  = np.asarray(weir1).flatten()
                _culvert_boundaries['height'] = np.asarray(weir2).flatten()
                _culvert_boundaries['subcritical_flow_coefficient'] = np.asarray(weir3).flatten()
                _culvert_boundaries['supercritical_flow_coefficient'] = np.asarray(weir4).flatten()
                _culvert_boundaries['cross_barrier_pipe_height'] = np.asarray(weir5).flatten()
                _culvert_boundaries['friction_factor'] = np.asarray(weir6).flatten()
                _culvert_boundaries['pipe_diameter'] = np.asarray(weir7).flatten()
                _culvert_boundaries['btype'] = btype
                culvert_boundaries.append(_culvert_boundaries)
    except:
        ocean_boundaries = None
        land_boundaries = None
        inner_boundaries = None
        inflow_boundaries = None
        outflow_boundaries = None
        weir_boundaries = None
        culvert_boundaries = None
        
    if not ocean_boundaries: ocean_boundaries = None
    if not land_boundaries: land_boundaries   = None
    if not inner_boundaries: inner_boundaries = None
    if not inflow_boundaries: inflow_boundaries = None
    if not outflow_boundaries: outflow_boundaries = None
    if not weir_boundaries: weir_boundaries = None
    if not culvert_boundaries: culvert_boundaries = None
    grid = dict()
    grid['x']                  = x
    grid['y']                  = y
    grid['values']             = z
    grid['elements']           = elements
    grid['nodeID']             = nodeID
    grid['elementID']          = elementID
    grid['ocean_boundaries']   = ocean_boundaries
    grid['land_boundaries']    = land_boundaries
    grid['inner_boundaries']   = inner_boundaries
    grid['weir_boundaries']    = weir_boundaries
    grid['inflow_boundaries']  = inflow_boundaries
    grid['outflow_boundaries'] = outflow_boundaries
    grid['culvert_boundaries'] = culvert_boundaries
    grid['description']        = description
    f.close()
    return grid

  def make_plot(self, extent=None, epsg=None, axes=None, title=None, total_colors=256, cbar_label='elevation [m]', **kwargs):
    self._init_fig(axes, extent, title, epsg)
    vmin = kwargs.pop("vmin", np.min(self.values[self._idx]))
    vmax = kwargs.pop("vmax", np.max(self.values[self._idx]))
    show = kwargs.pop("show", False)
    if vmax < 0.:
      levels = kwargs.pop("levels", np.linspace(vmin, vmax, total_colors))
      mlevel = np.mean(self.values[self._idx])
      cmap   = kwargs.pop("cmap", plt.cm.seismic)
      col_val = 0.5      
    else:
      mlevel = 0.
      wet_count = int(np.floor(total_colors * \
              (float((self.values[self._idx] < 0.).sum())/self.values[self._idx].size)))
      col_val = float(wet_count)/float(total_colors)
      dry_count = total_colors - wet_count
      colors_undersea = plt.cm.bwr(np.linspace(1., 0., wet_count))
      colors_land = plt.cm.terrain(np.linspace(0.25, 1., dry_count))
      colors = np.vstack((colors_undersea, colors_land))
      cmap = LinearSegmentedColormap.from_list('cut_terrain', colors)
      wlevels = np.linspace(vmin, mlevel, wet_count, endpoint=False)
      dlevels = np.linspace(mlevel, vmax, dry_count)
      levels = np.hstack((wlevels, dlevels))
    norm = FixPointNormalize(sealevel=mlevel, vmax=vmax, vmin=vmin, col_val=col_val)
    self._axes.tricontourf(self.x, self.y, self.elements, self.values,
                           levels=levels, cmap=cmap, norm=norm, extend='both', **kwargs)
    self._init_cbar(cmap, vmin, vmax)
    self._cbar.set_label(cbar_label)
    self._cbar.set_ticks([vmin, vmin + col_val *(vmax-vmin), vmax])
    self._cbar.set_ticklabels([np.around(vmin, 2), mlevel, np.around(vmax, 2)])
    if show == True:
      plt.show()

  def interpolate_DEM(self, tile, **kwargs):
    method = kwargs.pop("method", "FVM")
    channel_polygons = kwargs.pop("channel_polygons", None)
    bar_polygons = kwargs.pop("bar_polygons", None)
    tile_extent = tile.get_extent()
    idxs_of_points_in_tile = self.get_extent_idx(tile_extent, tile.epsg)
    xyz = tile.get_xyz(epsg=self.epsg)
    for idx in idxs_of_points_in_tile:
      path = self.get_finite_volume_Path(idx)
      _idx, = np.where(np.logical_and(
                  np.logical_and(xyz[:,0]>=np.min(path.vertices[:,0]), xyz[:,0]<=np.max(path.vertices[:,0])),
                  np.logical_and(xyz[:,1]>=np.min(path.vertices[:,1]), xyz[:,1]<=np.max(path.vertices[:,1]))))
      values = xyz[_idx,:]
      if method=="FVM":
        _values = values[np.where(path.contains_points(values[:,0:2])),2]
        self.values[idx] = np.mean(_values)
      else:
        self.values[idx] = griddata((xyz[_idx,0], xyz[_idx,1]), xyz[_idx,2], (self.x[idx], self.y[idx]), method=method)
      if channel_polygons is not None:
        for channel_polygon in channel_polygons:
          if channel_polygon.contains_point((self.x[idx], self.y[idx])):#path.intersects_path(channel_polygon):
            try: 
              _values
            except:
              _values = values[np.where(path.contains_points(values[:,0:2])),2]
            self.values[idx] = np.min(_values)
  
      if bar_polygons is not None:
        for bar_polygon in bar_polygons:
          if bar_polygon.contains_point((self.x[idx], self.y[idx])):#path.intersects_path(bar):
            try:
              _values
            except:
              _values = values[np.where(path.contains_points(values[:,0:2])),2]
            self.values[idx] = np.max(_values)

  def write_fort14(self, path):
    if self.datum != 'MSL':
      self.values = self.values + self.datum_offset
    f = open(path, 'w')
    f.write(self.description + '\n')
    f.write("{}  {}\n".format(self.elements.shape[0], self.values.size))
    for i in range(self.values.size):
      f.write("{:>10} {:15.10f} {:15.10f} {:15.10f}\n".format(
            self.nodeID[i]+1, self.x[i], self.y[i], -self.values[i]))
    for i in range(self.elements.shape[0]):
      f.write("{:5d}{:5d} {:5d} {:5d} {:5d}\n".format(self.elementID[i], 3,self.elements[i,0]+1, self.elements[i,1]+1, self.elements[i,2]+1))
    if self.ocean_boundaries is None:
      ocean_boundaries = []
    else:
      ocean_boundaries = self.ocean_boundaries
    _sum=0
    for i in range(len(ocean_boundaries)):
      _sum +=  ocean_boundaries[i].size
    f.write("{:d} = Number of open boundaries\n".format(len(ocean_boundaries)))
    f.write("{:d} = Total number of open boundary nodes\n".format(_sum))
    if self.ocean_boundaries is not None:
      for i in range(len(self.ocean_boundaries)):
        f.write("{:d} = Number of nodes for open boundary {:d}\n".format(ocean_boundaries[i].size, i+1))
        for j in range(ocean_boundaries[i].size):
          f.write("{:d}\n".format(ocean_boundaries[i][j]+1))
    remainingBoundaries = 0
    _sum=0
    if self.land_boundaries is not None:
      remainingBoundaries+=len(self.land_boundaries)
      for i in range(len(self.land_boundaries)):
        _sum +=  self.land_boundaries[i][0].size
   
    if self.inner_boundaries is not None:
      remainingBoundaries+=len(self.inner_boundaries)
      for i in range(len(self.inner_boundaries)):
        _sum +=  self.inner_boundaries[i][0].size
    
    if self.inflow_boundaries is not None:
      remainingBoundaries+=len(self.inflow_boundaries)
      for i in range(len(self.inflow_boundaries)):
        _sum +=  self.inflow_boundaries[i][0].size
    
    if self.outflow_boundaries is not None:
      remainingBoundaries+=len(self.outflow_boundaries)
      for i in range(len(self.outflow_boundaries)):
        _sum +=  self.outflow_boundaries[i][0].size
            
    if self.weir_boundaries is not None:
      remainingBoundaries+=len(self.weir_boundaries)
      for i in range(len(self.weir_boundaries)):
        _sum +=  self.weir_boundaries[i]['front_face'].size
        _sum +=  self.weir_boundaries[i]['back_face'].size
    
    if self.culvert_boundaries is not None:
      remainingBoundaries+=len(self.culvert_boundaries)
      for i in range(len(self.culvert_boundaries)):
        _sum +=  self.culvert_boundaries[i]['front_face'].size
        _sum +=  self.culvert_boundaries[i]['back_face'].size
    
    f.write("{:d} = Number of land boundaries\n".format(remainingBoundaries))
    f.write("{:d} = Total number of land boundary nodes\n".format(_sum))
    
    _cnt = 1
    if self.land_boundaries is not None:
      for i in range(len(self.land_boundaries)):
        f.write("{:d} {:d} = Number of nodes for land boundary {:d}\n".format(
                self.land_boundaries[i][0].size, self.land_boundaries[i][-1], _cnt))
        _cnt+=1
        for j in range(self.land_boundaries[i][0].size):
          f.write("{:d}\n".format(self.land_boundaries[i][0][j]+1))
    
    if self.inner_boundaries is not None:
      for i in range(len(self.inner_boundaries)):
        f.write("{:d} {:d} = Number of nodes for closed (\"island\")  boundary (land boundary {:d})\n".format(
                    self.inner_boundaries[i][0].size, self.inner_boundaries[i][-1], _cnt))
        _cnt+=1
        for j in range(self.inner_boundaries[i][0].size):
          f.write("{:d}\n".format(self.inner_boundaries[i][0][j]+1))
        
    if self.inflow_boundaries is not None:
      for i in range(len(self.inflow_boundaries)):
        f.write("{:d} {:d} = Number of nodes for inflow boundary (land boundary {:d})\n".format(
                    self.inflow_boundaries[i][0].size, self.inflow_boundaries[i][-1], _cnt))
        _cnt+=1
        for j in range(self.inflow_boundaries[i][0].size):
          f.write("{:d}\n".format(self.inflow_boundaries[i][0][j]+1))
              
    if self.outflow_boundaries is not None:
      for i in range(len(self.outflow_boundaries)):
        f.write("{:d} {:d} = Number of nodes for outflow boundary (land boundary {:d})\n".format(
                    self.outflow_boundaries[i][0].size, self.outflow_boundaries[i][-1], _cnt))
        _cnt+=1
        for j in range(self.outflow_boundaries[i][0].size):
          f.write("{:d} {:.3f} {:.3f}\n".format(self.outflow_boundaries[i][0][j]+1, self.outflow_boundaries[i][1][j],
          self.outflow_boundaries[i][2][j]))
                
    if self.weir_boundaries is not None:
      for i in range(len(self.weir_boundaries)):
        f.write("{:d} {:d} = Number of node pairs for weir (land boundary {:d})\n".format(
                    self.weir_boundaries[i]['front_face'].size, self.weir_boundaries[i]['btype'], _cnt))
        _cnt+=1
        for j in range(self.weir_boundaries[i]['front_face'].size):
          f.write("{:d} {:d} {:.3f} {:.3f} {:.3f}\n".format(
          self.weir_boundaries[i]['front_face'][j]+1, self.weir_boundaries[i]['back_face'][j]+1,
          self.weir_boundaries[i]['height'][j], self.weir_boundaries[i]['subcritical_flow_coefficient'][j],
          self.weir_boundaries[i]['supercritical_flow_coefficient'][j]))
                
    if self.culvert_boundaries is not None:
      for i in range(len(self.culvert_boundaries)):
        f.write("{:d} {:d} = Number of nodes pairs for culvert boundary (land boundary {:d})\n".format(
                    len(self.culvert_boundaries[i]['front_face']), self.culvert_boundaries[i]['btype'], _cnt))
        _cnt+=1
        for j in range(self.culvert_boundaries[i]['front_face'].size):
          f.write("{:d} {:d} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f}\n".format(
          self.culvert_boundaries[i]['front_face'][j]+1, self.culvert_boundaries[i]['back_face'][j]+1,
          self.culvert_boundaries[i]['height'][j],     self.culvert_boundaries[i]['subcritical_flow_coefficient'][j],
          self.culvert_boundaries[i]['supercritical_flow_coefficient'][j],     self.culvert_boundaries[i]['cross_barrier_pipe_height'][j],
          self.culvert_boundaries[i]['friction_factor'][j],     self.culvert_boundaries[i]['pipe_diameter'][j]))
    f.close()
    if self.datum != 'MSL':
      self.values = self.values - self.datum_offset

  def _init_fort13(self, path, spinup_attributes, runtime_attributes):
    if path is None:
      self.fort13 = None
    else:
      self.fort13=NodalAttributes.parse_fort13(path, spinup_attributes, runtime_attributes)
