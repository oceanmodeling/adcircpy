from collections import OrderedDict
from AdcircPy.Tides import orbital_constants

def init_constituents(self):
  constituents = OrderedDict(sorted(orbital_constants.orbital_frequency.items(), key=lambda x: x[1]))
  TidalDB = OrderedDict()
  for constituent in constituents:
    TidalDB[constituent] = dict()
    TidalDB[constituent]['orbital_frequency'] = orbital_constants.orbital_frequency[constituent]
    if constituent in orbital_constants.doodson_coefficient.keys():
      TidalDB[constituent]['doodson_coefficient'] = orbital_constants.doodson_coefficient[constituent]
    if constituent in orbital_constants.tidal_potential_amplitude.keys():
      TidalDB[constituent]['tidal_potential_amplitude'] = orbital_constants.tidal_potential_amplitude[constituent]
    if constituent in orbital_constants.earth_tidal_potential_reduction_factor.keys():
      TidalDB[constituent]['earth_tidal_potential_reduction_factor'] = orbital_constants.earth_tidal_potential_reduction_factor[constituent]
  for key in TidalDB:
    self[key] = TidalDB[key]
