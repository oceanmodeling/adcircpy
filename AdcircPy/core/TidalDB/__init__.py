import os
import numpy as np
import netCDF4
from AdcircPy.core.TidalDB._TidalDB import orbital_frequency,\
                             doodson_coefficient,\
                             tidal_potential_amplitude,\
                             earth_tidal_potential_reduction_factor

_TPXO_default_path = os.path.dirname(os.path.abspath(__file__)) + '/h_tpxo9.v1.nc'

class TidalDB(object):
    def __init__(self):
        self.orbital_frequency                      = orbital_frequency
        self.doodson_coefficient                    = doodson_coefficient
        self.tidal_potential_amplitude              = tidal_potential_amplitude
        self.earth_tidal_potential_reduction_factor = earth_tidal_potential_reduction_factor

    def generate_equilibrium_factors(self, start_date, end_date):
        pass




class TPXO(object):
    def __init__(self):
        self.Dataset = netCDF4.Dataset(_TPXO_default_path)
        self.constituents = self.Dataset['con'][:].tostring().decode('UTF-8').split()

    def get_constituents_at_lonlat(self, qlon, qlat, constituent_list):
        return _TPXO.get_constituents_at_lonlat(self, qlon, qlat, constituent_list)


