import numpy as np
import utm
from haversine import haversine
from scipy.interpolate import interp1d
import pyproj
import matplotlib.pyplot as plt    
from AdcircPy.Winds._ATCF import _ATCF

class _WindVortex(_ATCF):

  def __init__(self, storm_id, start_date=None, end_date=None):
    super(_WindVortex, self).__init__(storm_id, start_date, end_date)
    self.__fix_background_pressure()
    self.__fix_radius_of_last_closed_isobar()
    self.__fix_radius_of_maximum_winds()
    self.__init_vars()
    self.__init_record_number()
    self.__init_speed_and_direction()
    self.__init_fort22()

  @property
  def u(self):
    return self._u

  @property
  def v(self):
    return self._v

  def printf(self):
    print('\n', end="")
    for line in self._fort22:
      print(line, end="")

  def plot_track(self):
    self.plt = plt
    self.plt.plot(self.longitude, self.latitude)
    self.plt.quiver(self.longitude, self.latitude, self.u, self.v)

  def __fill_missing_values(self, ATCF_list):
    if np.any(np.asarray(ATCF_list)==None):
      item = np.asarray(ATCF_list)
      item_x = np.arange(item.size)
      idx = np.where(item!=None)
      f = interp1d(item_x[idx], item[idx], fill_value="extrapolate")
      return [int(np.around(x,0)) for x in f(np.arange(item_x.size))]
    else:
      return ATCF_list

  def __fix_background_pressure(self):
    self._background_pressure = self.__fill_missing_values(self._background_pressure)

  def __fix_radius_of_last_closed_isobar(self):
    self._radius_of_last_closed_isobar = self.__fill_missing_values(self._radius_of_last_closed_isobar)

  def __fix_radius_of_maximum_winds(self):
    self._radius_of_maximum_winds = self.__fill_missing_values(self._radius_of_maximum_winds)

  def __init_vars(self):
    self._record_number = np.empty(len(self._datetime))
    self._direction = list()
    self._speed = list()
    self._number_of_isotachs = list()
    self._use_NEQ = list()
    self._use_SEQ = list()
    self._use_SWQ = list()
    self._use_NWQ = list()
    self._Rmax_of_NEQ = list()
    self._Rmax_of_SEQ = list()
    self._Rmax_of_SWQ = list()
    self._Rmax_of_NWQ = list()
    self._HollandB = list()
    self._HollandB_of_NEQ = list()
    self._HollandB_of_SEQ = list()
    self._HollandB_of_SWQ = list()
    self._HollandB_of_NWQ = list()
    self._Vmax_of_NEQ = list()
    self._Vmax_of_SEQ = list()
    self._Vmax_of_SWQ = list()
    self._Vmax_of_NWQ = list()
    self.__unique_datetimes = list(sorted(set(self._datetime)))
  
  def __init_record_number(self):
    datetime_dict = dict()
    for i, datetime in enumerate(self.__unique_datetimes):
      indexes = np.where(np.asarray(self._datetime)==datetime)
      for idx in indexes:
        self._record_number[idx]=i+1
    self._record_number = [int(x) for x in self._record_number]

  def __init_uv(self):
    self._u = list()
    self._v = list()
    for speed, direction in zip(self._speed, self._direction):
      direction=np.deg2rad(np.abs((450-direction))%360)
      self._u.append(speed*np.cos(direction))
      self._v.append(speed*np.sin(direction))

  def __init_speed_and_direction(self):
    zone = utm.from_latlon(self._latitude[0], self._longitude[0])[2]
    utm_proj = pyproj.Proj(proj='utm', zone=zone)
    x, y = utm_proj(self._longitude, self._latitude)
    for i, datetime in enumerate(self.__unique_datetimes):
      indexes, = np.where(np.asarray(self.datetime)==datetime)
      for idx in indexes:
        if indexes[-1]+1<len(self.datetime):
          dt = (self.datetime[indexes[-1]+1] - self.datetime[idx]).total_seconds()/(60.*60.)
          dx = haversine((self._latitude[idx], self._longitude[indexes[-1]+1]), (self._latitude[idx], self._longitude[idx]),nautical_miles=True)
          dy = haversine((self._latitude[indexes[-1]+1], self._longitude[idx]), (self._latitude[idx], self._longitude[idx]),nautical_miles=True)
          vx = np.copysign(dx/dt, self._longitude[indexes[-1]+1]-self._longitude[idx])
          vy = np.copysign(dy/dt, self._latitude[indexes[-1]+1]-self._latitude[idx])
        else:
          dt = (self.datetime[idx]-self.datetime[indexes[0]-1]).total_seconds()/(60.*60.)
          dx = haversine((self._latitude[idx], self._longitude[indexes[0]-1]), (self._latitude[idx], self._longitude[idx]),nautical_miles=True)
          dy = haversine((self._latitude[indexes[0]-1], self._longitude[idx]), (self._latitude[idx], self._longitude[idx]),nautical_miles=True)
          vx = np.copysign(dx/dt, self._longitude[idx] - self._longitude[indexes[0]-1])
          vy = np.copysign(dy/dt, self._latitude[idx] - self._latitude[indexes[0]-1])
        speed = np.sqrt(dx**2+dy**2)/dt
        bearing = (360. + np.rad2deg(np.arctan2(vx, vy))) % 360
        self._speed.append(int(np.around(speed, 0)))
        self._direction.append(int(np.around(bearing,0)))
    self.__init_uv()

  def __init_fort22(self):
    self._fort22=list()
    for i in range(len(self._datetime)):
      string = "{:<2},".format(self._basin[i])
      string+= "{:>3},".format(self._storm_number[i])
      string+= "{:>11},".format(self._datetime[i].strftime('%Y%m%d%H'))
      string+= "{:3},".format("")
      string+= "{:>5},".format(self._record_type[i])
      string+= "{:>4},".format(int((self.datetime[i]-self.start_date).total_seconds()/3600))
      if self._latitude[i]>=0:
        string+= "{:>4}N,".format(int(self._latitude[i]/.1))
      else:
        string+= "{:>4}S,".format(int(self._latitude[i]/-.1))
      if self._longitude[i]>=0:
        string+= "{:>5}E,".format(int(self._longitude[i]/.1))
      else:
        string+= "{:>5}W,".format(int(self._longitude[i]/-.1))
      string+= "{:>4},".format(int(self._max_sustained_wind_speed[i]))
      string+= "{:>5},".format(int(self._central_pressure[i]))
      string+= "{:>3},".format(self._development_level[i])
      string+= "{:>4},".format(int(self._isotach[i]))
      string+= "{:>4},".format(self._quadrant[i])
      string+= "{:>5},".format(int(self._radius_for_NEQ[i]))
      string+= "{:>5},".format(int(self._radius_for_SEQ[i]))
      string+= "{:>5},".format(int(self._radius_for_SWQ[i]))
      string+= "{:>5},".format(int(self._radius_for_NWQ[i]))
      if self._background_pressure[i] <= self._central_pressure[i] and \
           1013>self._central_pressure[i]:
        string+= "{:>5},".format(1013)
      elif self._background_pressure[i] <= self._central_pressure[i] and \
           1013<=self._central_pressure[i]:
        string+= "{:>5},".format(int(self._central_pressure[i]+1))
      else:
        string+= "{:>5},".format(int(self._background_pressure[i]))
      string+= "{:>5},".format(int(self._radius_of_last_closed_isobar[i]))
      string+= "{:>4},".format(int(self._radius_of_maximum_winds[i]))
      string+= "{:>5},".format('') # gust
      string+= "{:>4},".format('') # eye
      string+= "{:>4},".format('') # subregion
      string+= "{:>4},".format('') # maxseas
      string+= "{:>4},".format('') # initials
      string+= "{:>3},".format(self._direction[i])
      string+= "{:>4},".format(self._speed[i])
      string+= "{:^12},".format(self._name[i])
      # from this point forwards it's all aswip
      string+= "{:>4},".format(self._record_number[i])
      # string+= "{:>5},".format(self._record_number[i])
      string+="\n"
      self._fort22.append(string)