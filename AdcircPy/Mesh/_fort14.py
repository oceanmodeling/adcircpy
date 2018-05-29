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
    grid['fort14_description'] = description
    grid['fort14_path']        = path
    # grid['levels'] = np.negative(np.flip([-40,-30,-20,-10,-5,0,2,4,6,8,10,15,20,30,40,50,100,200,300,400,500,550], axis=0))
    return grid

def write_fort14(self, path):
    if self.datum != 'MSL':
        self.values = self.values + self.datum_offset
    f = open(path, 'w')
    f.write(self.fort14_description + '\n')
    f.write("{}  {}\n".format(self.elements.shape[0], self.values.size))
    for i in range(self.values.size):
        f.write("{:>10} {:15.10f} {:15.10f} {:15.10f}\n".format(
            self.nodeID[i]+1, self.x[i], self.y[i], -self.values[i]))
    for i in range(self.elements.shape[0]):
        f.write("{:5d}{:5d} {:5d} {:5d} {:5d}\n".format(self.elementID[i], 3,self.elements[i,0]+1, self.elements[i,1]+1, self.elements[i,2]+1))
    if self.oceanBoundaries is None:
        oceanBoundaries = []
    else:
        oceanBoundaries = self.oceanBoundaries
    _sum=0
    for i in range(len(oceanBoundaries)):
        _sum +=  oceanBoundaries[i].size
    f.write("{:d} = Number of open boundaries\n".format(len(oceanBoundaries)))
    f.write("{:d} = Total number of open boundary nodes\n".format(_sum))
    if self.oceanBoundaries is not None:
        for i in range(len(self.oceanBoundaries)):
            f.write("{:d} = Number of nodes for open boundary {:d}\n".format(oceanBoundaries[i].size, i+1))
            for j in range(oceanBoundaries[i].size):
                f.write("{:d}\n".format(oceanBoundaries[i][j]+1))
    remainingBoundaries = 0
    _sum=0
    if self.landBoundaries is not None:
        remainingBoundaries+=len(self.landBoundaries)
        for i in range(len(self.landBoundaries)):
            _sum +=  self.landBoundaries[i][0].size
   
    if self.innerBoundaries is not None:
        remainingBoundaries+=len(self.innerBoundaries)
        for i in range(len(self.innerBoundaries)):
            _sum +=  self.innerBoundaries[i][0].size
    
    if self.inflowBoundaries is not None:
        remainingBoundaries+=len(self.inflowBoundaries)
        for i in range(len(self.inflowBoundaries)):
            _sum +=  self.inflowBoundaries[i][0].size
    
    if self.outflowBoundaries is not None:
        remainingBoundaries+=len(self.outflowBoundaries)
        for i in range(len(self.outflowBoundaries)):
            _sum +=  self.outflowBoundaries[i][0].size
            
    if self.weirBoundaries is not None:
        remainingBoundaries+=len(self.weirBoundaries)
        for i in range(len(self.weirBoundaries)):
            _sum +=  self.weirBoundaries[i]['front_face'].size
            _sum +=  self.weirBoundaries[i]['back_face'].size
    
    if self.culvertBoundaries is not None:
        remainingBoundaries+=len(self.culvertBoundaries)
        for i in range(len(self.culvertBoundaries)):
            _sum +=  self.culvertBoundaries[i]['front_face'].size
            _sum +=  self.culvertBoundaries[i]['back_face'].size
    
    f.write("{:d} = Number of land boundaries\n".format(remainingBoundaries))
    f.write("{:d} = Total number of land boundary nodes\n".format(_sum))
    
    _cnt = 1
    if self.landBoundaries is not None:
        for i in range(len(self.landBoundaries)):
            f.write("{:d} {:d} = Number of nodes for land boundary {:d}\n".format(
                        self.landBoundaries[i][0].size, self.landBoundaries[i][-1], _cnt))
            _cnt+=1
            for j in range(self.landBoundaries[i][0].size):
                f.write("{:d}\n".format(self.landBoundaries[i][0][j]+1))
    
    if self.innerBoundaries is not None:
        for i in range(len(self.innerBoundaries)):
            f.write("{:d} {:d} = Number of nodes for closed (\"island\")  boundary (land boundary {:d})\n".format(
                        self.innerBoundaries[i][0].size, self.innerBoundaries[i][-1], _cnt))
            _cnt+=1
            for j in range(self.innerBoundaries[i][0].size):
                f.write("{:d}\n".format(self.innerBoundaries[i][0][j]+1))
        
    if self.inflowBoundaries is not None:
        for i in range(len(self.inflowBoundaries)):
            f.write("{:d} {:d} = Number of nodes for inflow boundary (land boundary {:d})\n".format(
                        self.inflowBoundaries[i][0].size, self.inflowBoundaries[i][-1], _cnt))
            _cnt+=1
            for j in range(self.inflowBoundaries[i][0].size):
                f.write("{:d}\n".format(self.inflowBoundaries[i][0][j]+1))
                
    if self.outflowBoundaries is not None:
        for i in range(len(self.outflowBoundaries)):
            f.write("{:d} {:d} = Number of nodes for outflow boundary (land boundary {:d})\n".format(
                        self.outflowBoundaries[i][0].size, self.outflowBoundaries[i][-1], _cnt))
            _cnt+=1
            for j in range(self.outflowBoundaries[i][0].size):
                f.write("{:d} {:.3f} {:.3f}\n".format(self.outflowBoundaries[i][0][j]+1, self.outflowBoundaries[i][1][j],
                self.outflowBoundaries[i][2][j]))
                
    if self.weirBoundaries is not None:
        for i in range(len(self.weirBoundaries)):
            f.write("{:d} {:d} = Number of node pairs for weir (land boundary {:d})\n".format(
                        self.weirBoundaries[i]['front_face'].size, self.weirBoundaries[i]['btype'], _cnt))
            _cnt+=1
            for j in range(self.weirBoundaries[i]['front_face'].size):
                f.write("{:d} {:d} {:.3f} {:.3f} {:.3f}\n".format(
                self.weirBoundaries[i]['front_face'][j]+1, self.weirBoundaries[i]['back_face'][j]+1,
                self.weirBoundaries[i]['height'][j], self.weirBoundaries[i]['subcritical_flow_coefficient'][j],
                self.weirBoundaries[i]['supercritical_flow_coefficient'][j]))
                
    if self.culvertBoundaries is not None:
        for i in range(len(self.culvertBoundaries)):
            f.write("{:d} {:d} = Number of nodes pairs for culvert boundary (land boundary {:d})\n".format(
                        len(self.culvertBoundaries[i]['front_face']), self.culvertBoundaries[i]['btype'], _cnt))
            _cnt+=1
            for j in range(self.culvertBoundaries[i]['front_face'].size):
                f.write("{:d} {:d} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f}\n".format(
                self.culvertBoundaries[i]['front_face'][j]+1, self.culvertBoundaries[i]['back_face'][j]+1,
                self.culvertBoundaries[i]['height'][j],     self.culvertBoundaries[i]['subcritical_flow_coefficient'][j],
                self.culvertBoundaries[i]['supercritical_flow_coefficient'][j],     self.culvertBoundaries[i]['cross_barrier_pipe_height'][j],
                self.culvertBoundaries[i]['friction_factor'][j],     self.culvertBoundaries[i]['pipe_diameter'][j]))
    f.close()

    if self.datum != 'MSL':
        self.values = self.values - self.datum_offset
