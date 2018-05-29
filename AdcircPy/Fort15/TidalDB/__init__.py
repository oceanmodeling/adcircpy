import os
import numpy as np
import netCDF4
from scipy.interpolate import griddata, RectBivariateSpline
# from scipy.io import FortranFile

TPXO_path = os.path.dirname(os.path.abspath(__file__)) + '/h_tpxo9.v1.nc'

class TidalDB():
    def __init__(self):
        # from wikipedia https://en.wikipedia.org/wiki/Theory_of_tides
        self.orbital_frequency = {  'M4'      : 0.0002810378050173,
                                    'M6'      : 0.0004215567080107,
                                    'MK3'     : 0.0002134400613513,
                                    'S4'      : 0.0002908882086657,
                                    'MN4'     : 0.0002783986019952,
                                    'S6'      : 0.0004363323129986,
                                    'M3'      : 0.0002107783537630,
                                    '2"MK3'   : 0.0002081166466594,
                                    'M8'      : 0.0005620756090649,
                                    'MS4'     : 0.0002859630068415,
                                    'M2'      : 0.0001405189025086,
                                    'S2'      : 0.0001454441043329,
                                    'N2'      : 0.0001378796994865,
                                    'Nu2'     : 0.0001382329037065,
                                    'MU2'     : 0.0001355937006844,
                                    '2"N2'    : 0.0001352404964644,
                                    'lambda2' : 0.0001428049013108,
                                    'T2'      : 0.0001452450073529,
                                    'R2'      : 0.0001456432013128,
                                    '2SM2'    : 0.0001503693061571,
                                    'L2'      : 0.0001431581055307,
                                    'K2'      : 0.0001458423172006,
                                    'K1'      : 0.0000729211583579,
                                    'O1'      : 0.0000675977441508,
                                    'OO1'     : 0.0000782445730498,
                                    'S1'      : 0.0000727220521664,
                                    'M1'      : 0.0000702594512543,
                                    'J1'      : 0.0000755603613800,
                                    'RHO'     : 0.0000653117453487,
                                    'Q1'      : 0.0000649585411287,
                                    '2Q1'     : 0.0000623193381066,
                                    'P1'      : 0.0000725229459750,
                                    'Mm'      : 0.0000026392030221,
                                    'Ssa'     : 0.0000003982128677,
                                    'Sa'      : 0.0000001991061914,
                                    'Msf'     : 0.0000049252018242,
                                    'Mf'      : 0.0000053234146919}
            
        self.doodson_coefficient = {'M4'     : 4,
                                        'M6'     : 6,
                                        'MK3'    : 3,
                                        'S4'     : 4,
                                        'MN4'    : 4,
                                        'S6'     : 6,
                                        'M3'     : 3,
                                        '2"MK3'  : 3,
                                        'M8'     : 8,
                                        'MS4'    : 4,
                                        'M2'     : 2,
                                        'S2'     : 2,
                                        'N2'     : 2,
                                        'Nu2'    : 2,
                                        'MU2'    : 2,
                                        '2"N2'   : 2,
                                        'lambda2': 2,
                                        'T2'     : 2,
                                        'R2'     : 2,
                                        '2SM2'   : 2,
                                        'L2'     : 2,
                                        'K2'     : 2,
                                        'K1'     : 1,
                                        'O1'     : 1,
                                        'OO1'    : 1,
                                        'S1'     : 1,
                                        'M1'     : 1,
                                        'J1'     : 1,
                                        'RHO'    : 1,
                                        'Q1'     : 1,
                                        '2Q1'    : 1,
                                        'P1'     : 1,
                                        'Mm'     : 0,
                                        'Ssa'    : 0,
                                        'Sa'     : 0,
                                        'Msf'    : 0,
                                        'Mf'     : 0}
        self.tidal_potential_amplitude = {
                                        'M2'  : 0.0001405189027,
                                        'S2'  : 0.0001454441043,  
                                        'N2'  : 0.0001378797074,  
                                        'K2'  : 0.0001458423017,  
                                        'K1'  : 0.141565,  
                                        'O1'  : 0.100514,  
                                        'P1'  : 0.0468,
                                        'Q1'  : 0.19256}


        self.earth_tidal_potential_reduction_factor = {
                                        'M2'  : 0.693,
                                        'S2'  : 0.693,  
                                        'N2'  : 0.693,  
                                        'K2'  : 0.693,  
                                        'K1'  : 0.736,  
                                        'O1'  : 0.695,  
                                        'P1'  : 0.706,
                                        'Q1'  : 0.695}

    def generate_equilibrium_factors(self, start_date, end_date):
        pass
        

class TPXO():
    def __init__(self):
        self.Dataset = netCDF4.Dataset(TPXO_path)
        self.constituents = self.Dataset['con'][:].tostring().decode('UTF-8').split()

    def get_constituents_at_lonlat(self, qlon, qlat, constituent_list):
        x = self.Dataset['lon_z'][:]
        y = self.Dataset['lat_z'][:]
        x = np.linspace(np.min(x),np.max(x),num=x.shape[0])
        y = np.linspace(np.min(y),np.max(y),num=y.shape[1])
        constituents = dict()
        for constituent in constituent_list:
            constituents[constituent] = dict()
            idx = self.constituents.index(constituent.lower())
            ha_interpolator = RectBivariateSpline(x, y, self.Dataset['ha'][idx,:,:])
            constituents[constituent]['ha'] = ha_interpolator.ev(qlon, qlat)
            hp_interpolator = RectBivariateSpline(x, y, self.Dataset['hp'][idx,:,:])
            constituents[constituent]['hp'] = hp_interpolator.ev(qlon, qlat)
        return constituents

