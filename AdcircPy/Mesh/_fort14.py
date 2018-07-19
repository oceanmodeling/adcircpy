import numpy as np

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
    # grid['fort14_path']        = path
    # grid['levels'] = np.negative(np.flip([-40,-30,-20,-10,-5,0,2,4,6,8,10,15,20,30,40,50,100,200,300,400,500,550], axis=0))
    return grid

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
