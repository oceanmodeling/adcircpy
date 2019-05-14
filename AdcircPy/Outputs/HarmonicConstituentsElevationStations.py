from collections import OrderedDict
from netCDF4 import Dataset
from AdcircPy.Tides import orbital_constants
from AdcircPy.Model._StationsOutput import _StationsOutput


class HarmonicConstituentsElevationStations(OrderedDict):

    def __init__(self, **station_data):
        super(HarmonicConstituentsElevationStations, self).__init__(
                                                                **station_data)

    @classmethod
    def from_ascii(cls, fort51, fort15):
        constituents = list()
        fort15 = _StationsOutput.parse_fort15(fort15, 'NOUTE')
        stations = dict()
        with open(fort51, 'r') as f:
            number_of_components = int(f.readline())
            for i in range(number_of_components):
                constituents.append(f.readline().split()[-1].strip())
            num_of_stations = int(f.readline())
            for station_id in range(num_of_stations):
                station_id = list(fort15.keys())[station_id]
                stations[station_id] = {'longitude': fort15[station_id]['x'],
                                        'latitude': fort15[station_id]['y'],
                                        'data': dict()}
                f.readline()
                for constituent in constituents:
                    line = f.readline().split()
                    amplitude = float(line[0])
                    phase = float(line[1])
                    stations[station_id]['data'][constituent] = {
                            "orbital_frequency":
                            orbital_constants.orbital_frequency[constituent],
                            "amplitude": amplitude,
                            "phase": phase
                            }
        return cls(**stations)

    @classmethod
    def from_netcdf(cls, fort51):
        nc = Dataset(fort51)
        stations = OrderedDict()
        for station in nc.variables['station_name'][:]:
            stations[station.tostring().strip().decode('utf-8')] = dict()
        for i, station in enumerate(stations.keys()):
            stations[station] = {'longitude': float(nc['x'][i]),
                                 'latitude': float(nc['y'][i]),
                                 'data': dict()}
            for j, component in enumerate(nc.variables['const'][:]):
                stations[station]['data'][component.tostring().strip().decode(
                            'utf-8')] = {
                                'orbital_frequency': float(nc['frequency'][j]),
                                'amplitude': float(nc['amp'][i, j].data),
                                'phase': float(nc['phs'][i, j].data)}
        return cls(**stations)
