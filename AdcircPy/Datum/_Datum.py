from __future__ import absolute_import, division, print_function
import numpy as np
from scipy.interpolate import griddata

            
def msl_to_navd88(self, datum_grid):
    if self.datum != 'MSL':
        return
    self.datum_grid=datum_grid
    f  = open(datum_grid) 
    f.readline().rstrip()
    NE, NP = map(int, f.readline().split())    
    lon    = list()
    lat    = list()
    offset  = list()

    for k in range(NP):
        line = f.readline().split()
        lon.append(float(line[1]))
        lat.append(float(line[2]))
        offset.append(float(line[3]))
    
    lon    = np.squeeze(lon)
    lat    = np.squeeze(lat)
    offset  = (-1)*np.squeeze(offset)

    if (np.all(lon == self.x)) \
            and (np.all(lat == self.y)) \
                and (offset.size == self.values.size):
        
        self.values = self.values - offset
        self.datum_offset = offset
    else:
        offset = griddata((lon, lat), offset, (self.x, self.y), fill_value=-99999.0)
        self.values = self.values - offset
        self.datum_offset = offset
    self.datum = 'NAVD88'

def navd88_to_msl(self):
    if self.datum == 'NAVD88':
        self.values = self.values + self.datum_offset
        delattr(self, 'datum_grid')
        delattr(self, 'datum_offset')
        self.datum = 'MSL'

