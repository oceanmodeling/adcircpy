import numpy as np
import copy
import warnings
from scipy.interpolate import griddata
from AdcircPy import Mesh as _Mesh
from AdcircPy import Datum
            
def msl_to_navd88(Datum_grid):
    f  = open(Datum_grid) 
    f.readline().rstrip()
    NE, NP = map(int, f.readline().split())    
    x    = list()
    y    = list()
    values  = list()
    for k in range(NP):
        line = f.readline().split()
        x.append(float(line[1]))
        y.append(float(line[2]))
        values.append(float(line[3]))
    x    = np.squeeze(x)
    y    = np.squeeze(y)
    values  = np.squeeze(values)
    f.close()
    kwargs =  { "x"      : x,
                "y"      : y,
                "values" : values,
                "_type"  : "navd88"}
    return Datum.DatumGrid(**kwargs)

def convert(self, Mesh, method='nearest'):

    if isinstance(Mesh, ("".__class__, u"".__class__)):
        Mesh = _Mesh.Mesh(Mesh)
    mesh = copy.deepcopy(Mesh)
    values = griddata((self.x, self.y), self.values, (Mesh.x, Mesh.y), method=method, fill_value=np.nan)
    Mesh.values = mesh.values + values
    Mesh.datum = self._type
    if np.any(np.isnan(Mesh.values)):
        warnings.warn("NaN values found during datum interpolation. Make sure the provided Datum conversion grid covers the entire domain.")
    return Mesh



