from collections import OrderedDict
from datetime import datetime, timedelta
import os
import calendar
import wget
import tarfile
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

def init_orbital_functions(self):
  self.orbital_functions_start = self._get_orbital_functions_start_of_record()
  self.orbital_functions_middle = self._get_orbital_functions_middle_of_record()
  
  print(self.orbital_functions_start["lunar_node_degrees"])
  print(self.orbital_functions_middle["lunar_node_degrees"])
  # self.mean_longitude_of_sun_start  = 280.1895014-.238724988*self.DYR+.9856473288*self.DDAY+.0410686387*self.spinup_date.hour
  # self.mean_longitude_of_sun_middle = 280.1895014-.238724988*self.DYR+.9856473288*self.DDAY+.0410686387*self.hour_middle
  # self.solar_perigee_start = 281.2208569+.01717836*self.DYR+.000047064*self.DDAY+.000001961*self.spinup_date.hour
  # self.solar_perigee_middle = 281.2208569+.01717836*self.DYR+.000047064*self.DDAY+.000001961*self.hour_middle
  # self.mean_longitude_of_moon_start=
  # self.mean_longitude_of_moon_middle=277.0256206+129.38482032*self.DYR+13.176396768*self.DDAY+.549016532*self.hour_middle

def _get_orbital_functions_start_of_record(self):
  return { "lunar_node_degrees"   : self._get_lunar_node_degrees(self.spinup_date.hour),
           "lunar_perigee"        : self._get_lunar_perigee_degrees(self.spinup_date.hour)
           "lunar_mean_longitude" : self._get_lunar_mean_longitude(self.spinup_date.hour)
  }

def _get_orbital_functions_middle_of_record(self):
  return { "lunar_node_degrees"   : self._get_lunar_node_degrees(self.hour_middle),
           "lunar_perigee"        : self._get_lunar_perigee_degrees(self.hour_middle),
           "lunar_mean_longitude" : self._get_lunar_mean_longitude(self.hour_middle) 
  }  

def _get_lunar_node_degrees(self, hours):
  # DN
  return (259.1560564-19.328185764*self.DYR-.0529539336*self.DDAY-.0022064139*hours) % 360.

def _get_lunar_perigee_degrees(self, hours):
  # DP
  return (334.3837214+40.66246584*self.DYR+.111404016*self.DDAY+.004641834*hours) % 360.

def _get_lunar_mean_longitude(self, hours):
  # DS
  return (277.0256206+129.38482032*self.DYR+13.176396768*self.DDAY+.549016532*hours) % 360.