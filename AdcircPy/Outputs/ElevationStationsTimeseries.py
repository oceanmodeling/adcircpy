from datetime import datetime, timedelta
from netCDF4 import Dataset
from AdcircPy.Outputs._OutputStations import _OutputStations


class ElevationStationsTimeseries(_OutputStations):

    def __init__(self, **stations):
        super(ElevationStationsTimeseries, self).__init__(**stations)

    @classmethod
    def from_netcdf(cls, path):
        nc = Dataset(path)
        time = nc['time'].base_date
        time = time.split('!')
        time = time[0]
        time = time.strip(' UTC')
        start_datetime = datetime.strptime(time, '%Y-%m-%d %H:%M:%S')
        time = [start_datetime + timedelta(seconds=x) for x in nc['time'][:]]
        station_ids = list()
        for station in nc.variables['station_name'][:]:
            station_ids.append(station.tostring().strip().decode('utf-8'))
        stations = dict()
        for i, station_id in enumerate(station_ids):
            stations[station_id] = {'x': nc['x'][i],
                                    'y': nc['y'][i],
                                    'values': nc['zeta'][:, i],
                                    'time': time}
        return cls(**stations)

    # def make_plot(self, station, **kwargs):
    #     axes = kwargs.pop('axes', None)
    #     title = kwargs.pop('title', None)
    #     show = kwargs.pop('show', False)
    #     legend = kwargs.pop('legend', True)

    #     if axes is None:
    #         fig = plt.figure()
    #         axes = fig.add_subplot(111)

    #     if title is not None:
    #         axes.set_title(title)

    #     idx = self.station_name.index(station)
    #     axes.plot(self.time, self.zeta[:, idx], **kwargs)

    #     if show is True:
    #         plt.legend()
    #         plt.show()
