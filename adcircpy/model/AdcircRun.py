from datetime import datetime, timedelta
from pathlib import Path
import numpy as np
from adcircpy.lib._Fort15 import _Fort15
from adcircpy.mesh import AdcircMesh as _AdcircMesh
from adcircpy.model import TidalForcing as _TidalForcing
from adcircpy.model.winds._WindForcing import _WindForcing


class AdcircRun(_Fort15):

    def __init__(
        self,
        start_date=None,
        end_date=None,
        netcdf=True,
        spinup_time=timedelta(0.)
    ):
        self.__start_date = start_date
        self.__end_date = end_date
        self.__spinup_time = timedelta(0.)
        self.__netcdf = netcdf
        self.__container = dict()
        self.__init_station_outputs()
        self.__init_global_outputs()
        super(AdcircRun, self).__init__()

    def set_elevation_stations_output(
        self,
        sampling_frequency,
        spinup=False,
        netcdf=True,
        harmonic_analysis=False
    ):
        self.__set_output_request(
            'stations', 'elevation', sampling_frequency, spinup, netcdf,
            harmonic_analysis)

    def set_velocity_stations_output(
        self,
        sampling_frequency,
        spinup=False,
        netcdf=True,
        harmonic_analysis=False
    ):
        self.__set_output_request(
            'stations', 'velocity', sampling_frequency, spinup, netcdf,
            harmonic_analysis)

    def set_meteorological_stations_output(
        self,
        sampling_frequency,
        spinup=False,
        netcdf=True,
        harmonic_analysis=False
    ):
        self.__set_output_request(
            'stations', 'meteorological', sampling_frequency, spinup, netcdf,
            harmonic_analysis)

    def set_concentration_stations_output(
        self,
        sampling_frequency,
        spinup=False,
        netcdf=True,
        harmonic_analysis=False
    ):
        self.__set_output_request(
            'stations', 'concentration', sampling_frequency, spinup, netcdf,
            harmonic_analysis)

    def set_elevation_global_output(
        self,
        sampling_frequency,
        spinup=False,
        netcdf=True,
        harmonic_analysis=False
    ):
        self.__set_output_request(
            'globals', 'elevation', sampling_frequency, spinup, netcdf,
            harmonic_analysis)

    def set_velocity_global_output(
        self,
        sampling_frequency,
        spinup=False,
        netcdf=True,
        harmonic_analysis=False
    ):
        self.__set_output_request(
            'globals', 'velocity', sampling_frequency, spinup, netcdf,
            harmonic_analysis)

    def set_meteorological_global_output(
        self,
        sampling_frequency,
        spinup=False,
        netcdf=True,
        harmonic_analysis=False
    ):
        self.__set_output_request(
            'globals', 'meteorological', sampling_frequency, spinup, netcdf,
            harmonic_analysis)

    def set_concentration_global_output(
        self,
        sampling_frequency,
        spinup=False,
        netcdf=True,
        harmonic_analysis=False
    ):
        self.__set_output_request(
            'globals', 'concentration', sampling_frequency, spinup, netcdf,
            harmonic_analysis)

    def add_elevation_output_station(self, station_name, vertices):
        self.__add_output_station('elevation', station_name, vertices)

    def add_velocity_output_station(self, station_name, vertices):
        self.__add_output_station('velocity', station_name, vertices)

    def add_meteorological_output_station(self, station_name, vertices):
        self.__add_output_station('meteorological', station_name, vertices)

    def add_concentration_output_station(self, station_name, vertices):
        self.__add_output_station('concentration', station_name, vertices)

    def remove_elevation_output_station(self, station_name):
        self.__remove_station('elevation', station_name)

    def remove_velocity_output_station(self, station_name):
        self.__remove_station('velocity', station_name)

    def remove_meteorological_output_station(self, station_name):
        self.__remove_station('meteorological', station_name)

    def remove_concentration_output_station(self, station_name):
        self.__remove_station('concentration', station_name)

    def get_elevation_global_output(self):
        return self.__get_output_request('globals', 'elevation')

    def get_velocity_global_output(self):
        return self.__get_output_request('globals', 'velocity')

    def get_meteorological_global_output(self):
        return self.__get_output_request('globals', 'meteorological')

    def get_concentration_global_output(self):
        return self.__get_output_request('globals', 'concentration')

    def get_elevation_output_stations(self):
        return self.__get_output_request('stations', 'elevation')

    def get_velocity_output_stations(self):
        return self.__get_output_request('stations', 'velocity')

    def get_meteorological_output_stations(self):
        return self.__get_output_request('stations', 'meteorological')

    def get_concentration_output_stations(self):
        return self.__get_output_request('stations', 'concentration')

    def dump(
        self,
        output_directory,
        overwrite=False,
        fort14='fort.14',
        fort13='fort.13',
        fort22='fort.22.best_track',
        coldstart='fort.15.coldstart',
        hotstart='fort.15.hotstart'
    ):
        output_directory = str(Path(output_directory))
        if fort14:
            path = str(Path(output_directory + '/' + fort14))
            self.mesh.write_fort14(path, overwrite)
        if fort13:
            path = str(Path(output_directory + '/' + fort13))
            self.mesh.write_fort13(path, overwrite)
        if fort22:
            if self.wind_forcing is not None:
                path = str(Path(output_directory + '/' + fort22))
                self.wind_forcing.dump(path, overwrite)
        if coldstart:
            path = str(Path(output_directory + '/' + coldstart))
            self.write_fort15('coldstart', path, overwrite)
        if hotstart:
            path = str(Path(output_directory + '/' + hotstart))
            self.write_fort15('hotstart', path, overwrite)

    def copy_fort15_stations(self, fort15):
        station_types = ['NOUTE', 'NOUTV', 'NOUTM', 'NOUTC']
        for station_type in station_types:
            stations = _Fort15.parse_stations(fort15, station_type)
            for name, vertices in stations.items():
                if station_type == 'NOUTE':
                    self.add_elevation_output_station(name, vertices)
                elif station_type == 'NOUTV':
                    self.add_velocity_output_station(name, vertices)
                elif station_type == 'NOUTM':
                    self.add_meteorological_output_station(name, vertices)
                elif station_type == 'NOUTC':
                    self.add_concentration_output_station(name, vertices)

    def write_fort15(self, runtype, path, overwrite=False):
        self.runtype = runtype
        if isinstance(path, bool) and path is True:
            print(self.fort15)
        else:
            if overwrite:
                with open(str(Path(path)), 'w') as f:
                    f.write(self.fort15)

    def __init_global_outputs(self):
        self.__container['globals'] = dict()
        self.__container['globals']['elevation'] = {
            'sampling_frequency': None,
            'spinup': False,
            'netcdf': True,
            'harmonic_analysis': False}
        self.__container['globals']['velocity'] = {
            'sampling_frequency': None,
            'spinup': False,
            'netcdf': True,
            'harmonic_analysis': False}
        self.__container['globals']['meteorological'] = {
            'sampling_frequency': None,
            'spinup': False,
            'netcdf': True,
            'harmonic_analysis': False}
        self.__container['globals']['concentration'] = {
            'sampling_frequency': None,
            'spinup': False,
            'netcdf': True,
            'harmonic_analysis': False}

    def __init_station_outputs(self):
        self.__container['stations'] = dict()
        self.__container['stations']['elevation'] = {
            'sampling_frequency': None,
            'spinup': False,
            'netcdf': True,
            'harmonic_analysis': False,
            'station_id': dict()}
        self.__container['stations']['velocity'] = {
            'sampling_frequency': None,
            'spinup': False,
            'netcdf': True,
            'harmonic_analysis': False,
            'station_id': dict()}
        self.__container['stations']['meteorological'] = {
            'sampling_frequency': None,
            'spinup': False,
            'netcdf': True,
            'harmonic_analysis': False,
            'station_id': dict()}
        self.__container['stations']['concentration'] = {
            'sampling_frequency': None,
            'spinup': False,
            'netcdf': True,
            'harmonic_analysis': False,
            'station_id': dict()}

    def __add_output_station(self, station_type, station_name, vertices):
        msg = "station_name {} ".format(station_name)
        msg += "already exists. Station names must be unique."
        msg += "{}".format(self.__container[
            'stations'][station_type]['station_id'].keys())
        assert station_name not in self.__container[
            'stations'][station_type]['station_id'].keys(), msg
        # assert self.mesh.has_station(vertices)
        self.__container['stations'][station_type][
            'station_id'][station_name] = vertices

    def __remove_station(self, station_type, station_name):
        self.__container['stations'][station_type][
            'station_id'].pop(station_name)

    def __set_output_request(
        self,
        output_type_1,
        output_type_2,
        sampling_frequency,
        spinup,
        netcdf,
        harmonic_analysis
    ):
        assert output_type_1 in self.__container.keys()
        assert output_type_2 in self.__container[output_type_1].keys()
        assert isinstance(sampling_frequency, (timedelta, type(None)))
        assert isinstance(spinup, bool)
        assert isinstance(netcdf, bool)
        assert isinstance(harmonic_analysis, bool)
        self.__container[output_type_1][output_type_2]['sampling_frequency'] \
            = sampling_frequency
        self.__container[output_type_1][output_type_2]['spinup'] = spinup
        self.__container[output_type_1][output_type_2]['netcdf'] = netcdf
        self.__container[output_type_1][output_type_2]['harmonic_analysis'] \
            = harmonic_analysis

    def __get_output_request(self, output_type_1, output_type_2):
        return self.__container[output_type_1][output_type_2]

    @property
    def mesh(self):
        try:
            return self.__mesh
        except AttributeError:
            raise AttributeError('Must set mesh attribute.')

    @property
    def adcirc_mesh(self):
        """ alias for mesh attribute """
        return self.mesh

    @property
    def tidal_forcing(self):
        try:
            if self.__start_date is not None:
                self.__tidal_forcing.start_date = self.start_date
            if self.__end_date is not None:
                self.__tidal_forcing.end_date = self.end_date
            self.__tidal_forcing.spinup_time = self.spinup_time
            return self.__tidal_forcing
        except AttributeError:
            return

    @property
    def wind_forcing(self):
        try:
            if self.__start_date is not None:
                self.__wind_forcing.start_date = self.start_date
            if self.__end_date is not None:
                self.__wind_forcing.end_date = self.end_date
            return self.__wind_forcing
        except AttributeError:
            return

    @property
    def start_date(self):
        if self.__start_date is not None:
            return self.__start_date
        else:
            raise AttributeError("Must set start_date of simulation.")

    @property
    def forcing_start_date(self):
        return self.start_date - self.spinup_time

    @property
    def spinup_time(self):
        try:
            return self.__spinup_time
        except AttributeError:
            if self.tidal_forcing is not None:
                self.tidal_forcing.spinup_time = timedelta(0.)
            return timedelta(0.)

    @property
    def spinup_factor(self):
        try:
            return self.__spinup_factor
        except AttributeError:
            return 1.

    @property
    def end_date(self):
        if self.__end_date is not None:
            if self.__start_date is not None:
                assert self.__end_date > self.__start_date
            return self.__end_date
        else:
            raise AttributeError("Must set end_date of simulation.")

    @property
    def netcdf(self):
        return self.__netcdf

    @start_date.setter
    def start_date(self, start_date):
        assert isinstance(start_date, (datetime, type(None)))
        if self.__end_date is not None:
            assert start_date < self.__end_date
        self.__start_date = start_date

    @end_date.setter
    def end_date(self, end_date):
        assert isinstance(end_date, (datetime, type(None)))
        if self.__start_date is not None:
            assert end_date > self.__start_date
        self.__end_date = end_date

    @spinup_time.setter
    def spinup_time(self, spinup_time):
        assert isinstance(spinup_time, timedelta)
        self.__spinup_time = np.abs(spinup_time)
        if self.tidal_forcing is not None:
            self.tidal_forcing.spinup_time = self.__spinup_time

    @spinup_factor.setter
    def spinup_factor(self, spinup_factor):
        spinup_factor = float(spinup_factor)
        assert (spinup_factor <= 0 and spinup_factor <= 1.)
        self.__spinup_factor = np.abs(spinup_factor)

    @netcdf.setter
    def netcdf(self, netcdf):
        assert isinstance(netcdf, bool)
        self.__netcdf = netcdf

    @mesh.setter
    def mesh(self, AdcircMesh):
        assert isinstance(AdcircMesh, _AdcircMesh)
        if AdcircMesh.SpatialReference is None:
            raise RuntimeError(
                'AdcircMesh._SpatialReference must be set before AdcircRun '
                + 'instantiation.')
        self.__mesh = AdcircMesh

    @adcirc_mesh.setter
    def adcirc_mesh(self, AdcircMesh):
        self.mesh = AdcircMesh

    @tidal_forcing.setter
    def tidal_forcing(self, TidalForcing):
        assert isinstance(TidalForcing, _TidalForcing)
        self.__tidal_forcing = TidalForcing

    @wind_forcing.setter
    def wind_forcing(self, WindForcing):
        assert isinstance(WindForcing, _WindForcing)
        self.__wind_forcing = WindForcing
