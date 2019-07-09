from collections import OrderedDict
from AdcircPy import Model
from AdcircPy import Outputs

def _from_ascii(path, fort14=None, fort15=None, datum=None, epsg=None, station_data=None):
    f = open(path)
    num_freq = int(f.readline())
    data_dict = OrderedDict()
    for i in range(num_freq):
        l = f.readline().split()
        data_dict[l[3]] = dict()
        data_dict[l[3]]['harmonic_frequency'] = float(l[0])
        data_dict[l[3]]['nodal_factor'] = float(l[1])
        data_dict[l[3]]['equilibrium_argument'] = float(l[2])

    num_datasets = int(f.readline())
    for i in range(num_datasets):
        _point = int(f.readline())
        for constituent in data_dict.keys():
            l = f.readline().split()
            data_dict[constituent][_point] = dict()
            data_dict[constituent][_point]['amplitude'] = float(l[0])
            data_dict[constituent][_point]['phase']     = float(l[1])
    
    if fort14==None and fort15==None and station_data == None:
        return HarmonicConstituents(data_dict)

    if fort14 is not None:
        if isinstance(fort14, Mesh.AdcircMesh):
            pass
        else:
            fort14 = Mesh.AdcircMesh.from_fort14(fort14)

    if station_data is not None:
        pass
    else:
        if fort15 is not None:
            if isinstance(fort15, Mesh.fort15):
                pass

