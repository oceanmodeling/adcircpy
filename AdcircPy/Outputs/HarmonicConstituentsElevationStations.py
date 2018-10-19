from AdcircPy.Tides import orbital_constants
from AdcircPy.Outputs._HarmonicConstituentsStations import _HarmonicConstituentsStations


class HarmonicConstituentsElevationStations(_HarmonicConstituentsStations):
  
  def __init__(self, fort51, fort15):
    super(HarmonicConstituentsElevationStations, self).__init__(fort15, 'NOUTE')
    self._fort51 = fort51
    self._init_fort51()

  def _init_fort51(self):
    self.constituents = list()
    with open(self._fort51, 'r') as f:
      number_of_components = int(f.readline())
      for i in range(number_of_components):
        self.constituents.append(f.readline().split()[-1].strip())
      num_of_stations = int(f.readline())
      for station_id in range(num_of_stations):
        station_id = list(self._fort15.keys())[station_id]
        self[station_id]={'longitude' : self._fort15[station_id]['x'],
                          'latitude'  : self._fort15[station_id]['y']}
        f.readline()
        for constituent in self.constituents:
          line = f.readline().split()
          amplitude = float(line[0])
          phase = float(line[1])
          self[station_id][constituent] = {"orbital_frequency"  : orbital_constants.orbital_frequency[constituent],
                                          "amplitude"         : amplitude,
                                          "phase"             : phase}
