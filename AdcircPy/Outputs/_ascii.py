from AdcircPy import Mesh
from AdcircPy.Outputs import _Outputs

def read_ascii_output(path, fort14=None, datum='MSL', epsg=4326, output_type=None):
    
    # Error checking for input args
    if fort14 is None:
        raise IOError("A fort.14 file is required to parse ASCII outputs.")
    if isinstance(fort14, ("".__class__, u"".__class__)):
        fort14 = Mesh.Mesh.init_from_fort14(fort14, datum, epsg)
    elif isinstance(fort14, Mesh.Mesh):
        pass
    else:
        raise IOError("fort14 keyword provided is neither a path to a fort.14 ASCII file nor an AdcircPy.Mesh instance!")

    f = open(path)
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
            if record_type == 1:
                _values.append(float(f.readline().split()[1]))
            elif record_type == 2:
                line = f.readline().split()
                _values.append((float(line[1]),float(line[2])))
            elif record_type == 3:
                line = f.readline().split()
                _values.append((float(line[1]), float(line[2]), float(line[3])))
        _values = np.asarray(_values)
        _values = np.ma.masked_equal(_values, -99999.0)
        values.append(_values)
    f.close()
    
    