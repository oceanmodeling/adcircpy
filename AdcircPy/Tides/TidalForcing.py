from collections import OrderedDict
from datetime import datetime, timedelta
import os
import calendar
import wget
import tarfile
import numpy as np
from netCDF4 import Dataset
from AdcircPy.Tides import orbital_constants

def get_cache_dir():
  cachedir = os.getenv('LOCALAPPDATA')
  if cachedir is None:
    cachedir = os.getenv('HOME')+'/.cache/AdcircPy'
  else: 
    cachedir += '/AdcircPy'
  return cachedir

def init_TPXO_cache(cachedir):
  os.makedirs(cachedir, exist_ok=True)
  if os.path.isfile(cachedir+"/h_tpxo9.v1.nc")==False:
    print('Building TPXO database cache on {}, please wait...'.format(cachedir+"/h_tpxo9.v1.nc"))
    print('(This will only happen the first time you run this software)')
    url='ftp://ftp.oce.orst.edu/dist/tides/Global/tpxo9_netcdf.tar.gz'
    if os.path.isfile(cachedir+"/tpxo9_netcdf.tar.gz")==False:
      wget.download(url, out=cachedir+"/tpxo9_netcdf.tar.gz")
    tpxo=tarfile.open(cachedir+"/tpxo9_netcdf.tar.gz")
    tpxo.extract('h_tpxo9.v1.nc', path=cachedir) # untar file into same directory
    tpxo.close()
    os.remove(cachedir+"/tpxo9_netcdf.tar.gz")

def init_constituent_dictionary(self):
  if self.constituents is not None:
    self.constituents = list(self.constituents)
    for constituent in self.constituents:
      if constituent not in list(self.keys()):
        raise Exception('\nUnknown Tidal Constituent \'{}\'.\n'.format(constituent)+\
                        'Possible constituents are: {} '.format(list(self.keys())))
  else:
    self.constituents = OrderedDict(sorted(orbital_constants.orbital_frequency.items(), key=lambda x: x[1]))
  for constituent in self.constituents:
    self[constituent] = dict()
    self[constituent]['orbital_frequency'] = orbital_constants.orbital_frequency[constituent]
    if constituent in orbital_constants.doodson_coefficient.keys():
      self[constituent]['doodson_coefficient'] = orbital_constants.doodson_coefficient[constituent]
    if constituent in orbital_constants.tidal_potential_amplitude.keys():
      self[constituent]['tidal_potential_amplitude'] = orbital_constants.tidal_potential_amplitude[constituent]
    if constituent in orbital_constants.earth_tidal_potential_reduction_factor.keys():
      self[constituent]['earth_tidal_potential_reduction_factor'] = orbital_constants.earth_tidal_potential_reduction_factor[constituent]

def init_spinup_date(self, spinup_date):
  if spinup_date is None:
    self.spinup_date = self.start_date - timedelta(days=15)
  elif isinstance(spinup_date, datetime):
    self.spinup_date = spinup_date
  else:
    raise IOError("spinup_date must be a datetime instance.")

def init_orbital_params(self):
  self.DYR  = self.spinup_date.year - 1900. 
  self.DDAY = self.spinup_date.timetuple().tm_yday + int((self.spinup_date.year-1901.)/4.)-1
  self.hour_middle = self.spinup_date.hour + ((self.end_date - self.spinup_date).total_seconds()/3600)/2
  self.N  = self.get_lunar_node(self.hour_middle)
  self.I  = np.arccos(.9136949-.0356926*np.cos(self.N))
  self.NU = np.arcsin(.0897056*np.sin(self.N)/np.sin(self.I))

def init_node_factors(self):
  for constituent in self.keys():
    # nodal factors are referenced to middle of record
    self[constituent]["nodal_factor"] = self._get_nodal_factor(constituent)
    # greenwich terms are referenced to the spinup_date
    self[constituent]["greenwich_term"] = self._get_greenwich_term(constituent)

def get_lunar_node(self, hours):
  # DN
  return np.deg2rad((259.1560564-19.328185764*self.DYR-.0529539336*self.DDAY-.0022064139*hours) % 360.)

def get_lunar_perigee(self, hours):
  # DP
  return np.deg2rad((334.3837214+40.66246584*self.DYR+.111404016*self.DDAY+.004641834*hours) % 360.)

def get_lunar_mean_longitude(self, hours):
  # DS
  return np.deg2rad((277.0256206+129.38482032*self.DYR+13.176396768*self.DDAY+.549016532*hours) % 360.)

def get_solar_perigee(self, hours):
  # DP1
  return np.deg2rad((281.2208569+.01717836*self.DYR+.000047064*self.DDAY+.000001961*hours) % 360.)

def get_solar_mean_longitude(self, hours):
  # DH
  return np.deg2rad((280.1895014-.238724988*self.DYR+.9856473288*self.DDAY+.0410686387*hours) % 360.)

def get_nodal_factor(self, constituent):
  if constituent   == "M2":
    return self._EQ78()
  elif constituent == "S2":
    return 1.0 # constant
  elif constituent == "N2":
    return self._EQ78()
  elif constituent == "K1":
    return self._EQ227()
  elif constituent == "M4":
    return (self._EQ78())**2.
  elif constituent == "O1":
    return self._EQ75()
  elif constituent == "M6":
    return (self._EQ78())**3.
  elif constituent == "MK3":
    return self._EQ78()*self._EQ227()
  elif constituent == "S4":
    return 1.0
  elif constituent == "MN4":
    return (self._EQ78())**2.

def _EQ78(self):
  return (np.cos(self.I/2)**4)/.91544
  
def _EQ227(self):
  return np.sqrt(.8965*np.sin(2.*self.I)**2+.6001*np.sin(2.*self.I)*np.cos(self.NU)+.1006)

def _EQ75(self):
  return np.sin(self.I)*np.cos(self.I/2.)**2/.37988

# def init_orbital_functions_start_of_record(self):
#     self.orbital_functions_start = dict()
#     DI =  
#     self.orbital_functions_start["I"] = 
#   # return { "lunar_node"           : self._get_lunar_node(self.spinup_date.hour),
#   #          "lunar_perigee"        : self._get_lunar_perigee(self.spinup_date.hour),
#   #          "lunar_mean_longitude" : self._get_lunar_mean_longitude(self.spinup_date.hour),
#   #          "solar_perigee"        : self._get_solar_perigee(self.spinup_date.hour),
#   #          "solar_mean_longitude" : self._get_solar_mean_longitude(self.spinup_date.hour)


# def init_orbital_functions_middle_of_record(self):
#   return { "lunar_node"           : self._get_lunar_node(self.hour_middle),
#            "lunar_perigee"        : self._get_lunar_perigee(self.hour_middle),
#            "lunar_mean_longitude" : self._get_lunar_mean_longitude(self.hour_middle),
#            "solar_perigee"        : self._get_solar_perigee(self.hour_middle),
#            "solar_mean_longitude" : self._get_solar_mean_longitude(self.hour_middle)
#   }  