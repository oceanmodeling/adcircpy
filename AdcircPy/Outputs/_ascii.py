def _parse_gridded(path, kwargs):
    
    f  = open(path)
    description       = f.readline().strip()
    line              = f.readline().split()
    num_of_datasets   = int(line[0])
    NP                = int(line[1])
    output_time       = float(line[2])
    output_interval   = float(line[3])
    record_type       = int(line[4]) # 1 for elevation, #2 for velocity #3 for 3D

    values   = list()
    time     = list()
    timestep = list()
    for dataset in range(num_of_datasets):
        line = f.readline().split()
        time.append(float(line[0]))
        timestep.append(int(line[1]))
        _values = list()
        for i in range(NP):
            _values.append(float(f.readline().split()[-1]))
        _values = np.asarray(_values).flatten()
        _values = np.ma.masked_equal(_values, -99999.0)
        values.append(_values)
    f.close()
    kwargs['values']   = values
    kwargs['time']     = time
    kwargs['timestep'] = timestep  
    return kwargs

def read_ascii_output(path, fort14=None, fort13=None, fort15=None, datum='MSL', epsg=4326, output_type=None):
    
    if fort14 is None:
        raise IOError("A fort.14 file is required to parse ASCII outputs.")

    f = open(path)
    header = f.readline()
    number_of_lines = int(f.readline().split()[1])
    f.close()
    
    if isinstance(fort14, ("".__class__, u"".__class__)):
        f = open(path)
        header = f.readline()
        number_of_nodes = int(f.readline().split()[1])
        f.close()
        string=True
    
    elif type(fort14)==type(adcpy.adcirc.topobathy()):
        number_of_nodes = fort14.values.shape[0]
    else:
        raise IOError("fort14 keyword provided is neither a path to a fort.14 file nor an adcirc topobathy instance!")

    if number_of_nodes == number_of_lines:
        if string == True:
            fort14 = adcpy.read_mesh(fort14)
        kwargs = {  "x"                  : fort14.x,
                    "y"                  : fort14.y,
                    "elements"           : fort14.elements,
                    "nodeID"             : fort14.nodeID,
                    "elementID"          : fort14.elementID, 
                    "oceanBoundaries"    : fort14.oceanBoundaries,
                    "landBoundaries"     : fort14.landBoundaries,
                    "innerBoundaries"    : fort14.innerBoundaries,
                    "weirBoundaries"     : fort14.weirBoundaries,
                    "inflowBoundaries"   : fort14.inflowBoundaries,
                    "outflowBoundaries"  : fort14.outflowBoundaries,
                    "culvertBoundaries"  : fort14.culvertBoundaries,
                    "fort14_path"        : fort14.fort14_path}
        return adcpy.adcirc.outputSurface(**adcpy.adcirc.outputs._parse_gridded(path, kwargs))
        
    else:
        print("This is a timeseries station output.")
        
        

    # else:
    #     print("I'm a string!")
    # if output_type is None:
    #     output_type = str(path[-2:])
    
    # if output_type in ['63'] and fort14 is None:
    #     raise IOError("A fort.14 file is required for reading ASCII surface type files.")
    
    # if output_type == '63':
    #     grid = adcpy.read_mesh(fort14, fort13, fort15, datum, epsg)
    #     return adcpy.adcirc.outputs._read_fort63(grid, path)

    # else:
    #     raise Exception("Output type not recognized.")
