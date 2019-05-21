import abc
from collections import OrderedDict
from collections.abc import Mapping
from datetime import timedelta


class _StationsOutput(Mapping, metaclass=abc.ABCMeta):

    def __init__(self, sampling_frequency=timedelta(0), netcdf=True,
                 spinup=False, harmonic_analysis=False):
        self._sampling_frequency = sampling_frequency
        self._netcdf = netcdf
        self._spinup = spinup
        self._harmonic_analysis = harmonic_analysis

    def __call__(self, UnstructuredMesh, runtype, forcing_start_date,
                 start_date, end_date, DTDP):
        return self.__get_fort15_params(
                        UnstructuredMesh, runtype, forcing_start_date,
                        start_date, end_date, DTDP)

    def __getitem__(self, key):
        return self._storage[key]

    def __iter__(self):
        return iter(self._storage)

    def __len__(self):
        return len(self._storage.keys())

    def add_station(self, x, y, station_id):
        self._storage[station_id] = (x, y)

    def _get_NHA(self, runtype):
        assert runtype in ['coldstart', 'hotstart', 'metonly']
        if runtype == 'coldstart':
            if self.spinup and self.harmonic_analysis:
                return 1
            else:
                return 0
        elif runtype == 'hotstart':
            if self.harmonic_analysis:
                return 1
            else:
                return 0
        else:
            NotImplementedError('Can we do harmonic analysis on winds only for'
                                + 'metonly runs?')
            return 0

    @abc.abstractmethod
    def add_stations_from_fort15(self, path, _hint):
        """
        Call this one using super() and pass _hint.
        """
        stations = self.__parse_stations_from_fort15(path, _hint)
        for station_id, (x, y) in stations.items():
            self._storage[station_id] = (x, y)

    def __get_fort15_params(self, UnstructuredMesh, runtype,
                            forcing_start_date, start_date, end_date, DTDP):
        NSPOOL = int(self.sampling_frequency.total_seconds()/DTDP)
        if runtype == 'coldstart':
            if self.spinup:
                if NSPOOL > 0.:
                    if self.netcdf:
                        NOUT = -5
                    else:
                        NOUT = -1
                    TOUTST = 0.
                    dt = start_date - forcing_start_date
                    TOUTF = dt.total_seconds() / (24.*60.*60.)
                    stations = self._storage
                else:
                    NOUT = 0
                    TOUTST = 0.
                    TOUTF = 0.
                    NSPOOL = 0
                    stations = dict()
            else:
                NOUT = 0
                TOUTST = 0.
                TOUTF = 0.
                NSPOOL = 0
                stations = dict()

        elif runtype == 'hotstart':
            if self.sampling_frequency.total_seconds() > 0.:
                if self.netcdf:
                    NOUT = -5
                else:
                    NOUT = -1
            else:
                NOUT = 0
            dt = start_date - forcing_start_date
            TOUTST = dt.total_seconds() / (24.*60.*60.)
            dt = end_date - forcing_start_date
            TOUTF = dt.total_seconds() / (24.*60.*60.)
            if NSPOOL > 0.:
                if self.netcdf:
                    NOUT = -5
                else:
                    NOUT = -1
            else:
                NOUT = 0
                TOUTST = 0.
                TOUTF = 0.
                NSPOOL = 0
                stations = dict()
            stations = self._storage
        pop_list = list()
        for station_id, (x, y) in stations.items():
            for InnerRing in UnstructuredMesh.InnerRings:
                if InnerRing.contains_point((x, y)):
                    pop_list.append(station_id)
                    continue
            if not UnstructuredMesh.OuterRing.contains_point((x, y)):
                pop_list.append(station_id)
        for station in pop_list:
            stations.pop(station)
        return NOUT, TOUTST, TOUTF, NSPOOL, stations

    @classmethod
    def _from_fort15(cls, path, _hint, sampling_frequency=timedelta(0),
                     netcdf=True, spinup=False, harmonic_analysis=False):
        cls = cls(sampling_frequency, netcdf, spinup, harmonic_analysis)
        stations = cls.__parse_stations_from_fort15(path, _hint)
        for station_id, (x, y) in stations.items():
            cls.add_station(x, y, station_id)
        return cls

    @staticmethod
    def __parse_stations_from_fort15(path, _hint):
        stations = OrderedDict()
        with open(path, 'r') as f:
            for line in f:
                if _hint in line:
                    line = f.readline().split('!')[0]
                    num = int(line)
                    for i in range(num):
                        line = f.readline().split('!')
                        if len(line) > 0:
                            station_name = line[1].strip(' \n')
                        else:
                            station_name = str(i)
                        coords = line[0].split(' ')
                        coords = [float(x) for x in coords if x != '']
                        x = float(coords[0])
                        y = float(coords[1])
                        stations[station_name] = (x, y)
        return stations

    @property
    @abc.abstractmethod
    def _storage(self):
        return self.__storage

    @property
    def stations(self):
        return self._storage.keys()

    @property
    def sampling_frequency(self):
        return self._sampling_frequency

    @property
    def spinup(self):
        return self._spinup

    @property
    def harmonic_analysis(self):
        return self._harmonic_analysis

    @property
    def netcdf(self):
        return self._netcdf

    @property
    def _sampling_frequency(self):
        return self.__sampling_frequency

    @property
    def _spinup(self):
        return self.__spinup

    @property
    def _harmonic_analysis(self):
        return self.__harmonic_analysis

    @property
    def _netcdf(self):
        return self.__netcdf

    @_sampling_frequency.setter
    def _sampling_frequency(self, sampling_frequency):
        assert isinstance(sampling_frequency, timedelta)
        self.__sampling_frequency = sampling_frequency

    @_spinup.setter
    def _spinup(self, spinup):
        assert isinstance(spinup, bool)
        self.__spinup = spinup

    @_harmonic_analysis.setter
    def _harmonic_analysis(self, harmonic_analysis):
        assert isinstance(harmonic_analysis, bool)
        self.__harmonic_analysis = harmonic_analysis

    @_netcdf.setter
    def _netcdf(self, netcdf):
        assert isinstance(netcdf, bool)
        self.__netcdf = netcdf
