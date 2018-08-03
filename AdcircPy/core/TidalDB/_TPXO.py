import numpy as np
from scipy.interpolate import RectBivariateSpline

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