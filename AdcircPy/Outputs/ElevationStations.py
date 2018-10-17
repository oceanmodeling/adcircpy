from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt 
from netCDF4 import Dataset
from AdcircPy.Outputs._OutputStations import _OutputStations


class ElevationStations(_OutputStations):
  def __init__(self, time, **station_data):
    super(ElevationStations, self).__init__(time, **station_data)
  
  @classmethod
  def from_netcdf(cls, path):
    nc = Dataset(path)
    try:
        start_datetime = datetime.strptime(nc['time'].units, 'seconds since %Y-%m-%d %H:%M:%S UTC')
        time = [start_datetime + timedelta(seconds=x) for x in nc['time'][:]]
    except:
        time = nc['time'][:]
    station_ids = list()
    for station in nc.variables['station_name'][:]:
        station_ids.append(station.tostring().strip().decode('utf-8'))
    station_data = dict()
    for i, station in enumerate(station_ids):
        station_data[station] = {
        'longitude' : nc['x'][i],
        'latitude'  : nc['y'][i],
        'zeta'      : nc['zeta'][:,i]}                                   
    return cls(time, **station_data)
    
  def make_plot(self, station, **kwargs):
    axes    = kwargs.pop('axes', None)
    title   = kwargs.pop('title', None)
    show    = kwargs.pop('show', False)
    legend  = kwargs.pop('legend', True)
    
    if axes is None:                
        fig = plt.figure()
        axes  = fig.add_subplot(111)

    if title is not None:
        axes.set_title(title)
    
    idx = self.station_name.index(station)
    axes.plot(self.time, self.zeta[:,idx], **kwargs)

    if show == True:
        plt.legend()
        plt.show()

