from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt 
from netCDF4 import Dataset
from AdcircPy import Outputs

def _from_netcdf(path):
    _Dataset      = Dataset(path)
    start_datetime = datetime.strptime(_Dataset['time'].units, 'seconds since %Y-%m-%d %H:%M:%S UTC')
    time = [start_datetime + timedelta(seconds=x) for x in _Dataset['time'][:]]
    station_ids = list()
    for station in _Dataset.variables['station_name'][:]:
        station_ids.append(station.tostring().strip().decode('utf-8'))
    station_data = dict()
    for i, station in enumerate(station_ids):
        station_data[station] = {
        'longitude' : _Dataset['x'][i],
        'latitude'  : _Dataset['y'][i],
        'zeta'      : _Dataset['zeta'][:,i]}                                   
    return Outputs.ElevationStations(time, **station_data)
    
def _make_plot(self, station, **kwargs):
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

