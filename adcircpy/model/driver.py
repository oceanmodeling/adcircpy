from datetime import datetime, timedelta
import subprocess
import pathlib
import shutil
import os
import numpy as np
from psutil import cpu_count
from adcircpy.model.fort15 import Fort15
from adcircpy.mesh.adcirc_mesh import AdcircMesh
from adcircpy.model.tidal_forcing import TidalForcing
from adcircpy.model.winds.wind_forcing import WindForcing


class AdcircRun(Fort15):

    def __init__(
        self,
        adcirc_mesh,
        tidal_forcing=None,
        wind_forcing=None,
        # wave_coupling=None,
        start_date=None,
        end_date=None,
        netcdf=True,
        spinup_time=None
    ):
        self._mesh = adcirc_mesh
        self._tidal_forcing = tidal_forcing
        self._wind_forcing = wind_forcing
        # self._wave_coupling = wave_coupling
        self._start_date = start_date
        self._end_date = end_date
        self._spinup_time = spinup_time
        self._netcdf = netcdf

    def add_elevation_output_station(self, station_name, vertices):
        self._elevation_stations[station_name] = self._certify_station(
            station_name, vertices, 'elevation')

    def add_velocity_output_station(self, station_name, vertices):
        self._velocity_stations[station_name] = self._certify_station(
            station_name, vertices, 'velocity')

    def add_meteorological_output_station(self, station_name, vertices):
        self._meteorological_stations[station_name] = self._certify_station(
            station_name, vertices, 'meteorological')

    def add_concentration_output_station(self, station_name, vertices):
        self._concentration_stations[station_name] = self._certify_station(
            station_name, vertices, 'concentration')

    def remove_elevation_output_station(self, station_name):
        self._elevation_stations.pop(station_name)

    def remove_velocity_output_station(self, station_name):
        self._velocity_stations.pop(station_name)

    def remove_meteorological_output_station(self, station_name):
        self._meteorological_stations.pop(station_name)

    def remove_concentration_output_station(self, station_name):
        self._concentration_stations.pop(station_name)

    def set_elevation_stations_output(
        self,
        sampling_frequency,
        spinup=False,
        netcdf=True,
        harmonic_analysis=False
    ):
        self._elevation_stations_output.update(
            self._certify_schema(
                sampling_frequency, spinup, netcdf, harmonic_analysis))

    def set_velocity_stations_output(
        self,
        sampling_frequency,
        spinup=False,
        netcdf=True,
        harmonic_analysis=False
    ):
        self._velocity_stations_output.update(
            self._certify_schema(
                sampling_frequency, spinup, netcdf, harmonic_analysis))

    def set_meteorological_stations_output(
        self,
        sampling_frequency,
        spinup=False,
        netcdf=True,
        harmonic_analysis=False
    ):
        self._meterological_stations_output.update(
            self._certify_schema(
                sampling_frequency, spinup, netcdf, harmonic_analysis))

    def set_concentration_stations_output(
        self,
        sampling_frequency,
        spinup=False,
        netcdf=True,
        harmonic_analysis=False
    ):
        self._concentration_stations_output.update(
            self._certify_schema(
                sampling_frequency, spinup, netcdf, harmonic_analysis))

    def set_elevation_surface_output(
        self,
        sampling_frequency,
        spinup=False,
        netcdf=True,
        harmonic_analysis=False
    ):
        self._elevation_surface_output.update(
            self._certify_schema(
                sampling_frequency, spinup, netcdf, harmonic_analysis))

    def set_velocity_surface_output(
        self,
        sampling_frequency,
        spinup=False,
        netcdf=True,
        harmonic_analysis=False
    ):
        self._velocity_surface_output.update(
            self._certify_schema(
                sampling_frequency, spinup, netcdf, harmonic_analysis))

    def set_meteorological_surface_output(
        self,
        sampling_frequency,
        spinup=False,
        netcdf=True,
        harmonic_analysis=False
    ):
        self._meteorological_surface_output.update(
            self._certify_schema(
                sampling_frequency, spinup, netcdf, harmonic_analysis))

    def set_concentration_surface_output(
        self,
        sampling_frequency,
        spinup=False,
        netcdf=True,
        harmonic_analysis=False
    ):
        self._concentration_surface_output.update(
            self._certify_schema(
                sampling_frequency, spinup, netcdf, harmonic_analysis))

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
        output_directory = pathlib.Path(output_directory).absolute()
        if fort14:
            path = output_directory / fort14
            self.mesh.write_fort14(path, overwrite)
        if fort13:
            path = output_directory / fort13
            self.mesh.write_fort13(path, overwrite)
        if fort22:
            if self.wind_forcing is not None:
                path = output_directory / fort22
                self.wind_forcing.dump(path, overwrite)
        if coldstart:
            path = output_directory / coldstart
            self.write_fort15('coldstart', path, overwrite)
        if hotstart:
            path = output_directory / hotstart
            self.write_fort15('hotstart', path, overwrite)

    def import_stations(self, fort15):
        station_types = ['NOUTE', 'NOUTV', 'NOUTM', 'NOUTC']
        for station_type in station_types:
            stations = Fort15.parse_stations(fort15, station_type)
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
        if isinstance(path, bool) and path is True:
            print(self.fort15(runtype))
        else:
            if overwrite:
                with open(str(pathlib.Path(path)), 'w') as f:
                    f.write(self.fort15(runtype))

    def adcirc(self, outdir, overwrite=False):
        outdir = pathlib.Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        self.dump(outdir, overwrite=overwrite)
        self._run_serial_coldstart(outdir, "adcirc")
        self._run_serial_hotstart(outdir, "adcirc")

    def padcirc(self, outdir, overwrite=False, nproc=-1, cmd='mpiexec'):
        assert cmd in ['mpiexec', 'srun']
        if cmd == 'mpiexec':
            nproc = cpu_count(logical=False) if nproc == -1 else nproc
            cmd += f" {nproc}"
            self._run_padcirc_coldstart(outdir, cmd)

    def _run_adcprep(self, rundir, nproc):
        wd = pathlib.Path.cwd()
        os.chdir(rundir.absolute())
        subprocess.check_call(['adcprep', '--np', f"{nproc:d}", '--partmesh'])
        subprocess.check_call(['adcprep', '--np', f"{nproc:d}", '--prepall'])
        os.chdir(wd.absolute())

    def _run_padcirc_coldstart(self, outdir, nproc, cmd):
        coldstart = outdir / 'coldstart'
        if coldstart.exists():
            shutil.rmtree(coldstart.absolute())
        coldstart.mkdir()
        os.symlink(
            outdir.absolute() / 'fort.14', coldstart.absolute() / 'fort.14')
        os.symlink(
            outdir.absolute() / 'fort.13', coldstart.absolute() / 'fort.13')
        os.symlink(
            outdir.absolute() / 'fort.15.coldstart',
            coldstart.absolute() / 'fort.15')
        self._run_adcprep(coldstart, nproc)
        wd = pathlib.Path.cwd()
        os.chdir(coldstart.absolute())
        subprocess.check_call([cmd, '-n', f"{nproc:d}", 'padcirc'])
        os.chdir(wd.absolute())

    def _run_serial_coldstart(self, outdir, cmd):
        coldstart = outdir / 'coldstart'
        if coldstart.exists():
            shutil.rmtree(coldstart.absolute())
        coldstart.mkdir()
        os.symlink(
            outdir.absolute() / 'fort.14', coldstart.absolute() / 'fort.14')
        os.symlink(
            outdir.absolute() / 'fort.13', coldstart.absolute() / 'fort.13')
        os.symlink(
            outdir.absolute() / 'fort.15.coldstart',
            coldstart.absolute() / 'fort.15')
        wd = pathlib.Path.cwd()
        os.chdir(coldstart.absolute())
        subprocess.check_call(cmd)
        os.chdir(wd.absolute())

    def _certify_station(self, station_name, vertices, physical_var):
        keys = list(self._container['stations'][physical_var].keys())
        vertices = np.asarray(vertices)
        msg = "vertices argument must be a two-tuple of floats."
        assert vertices.shape == (2,), msg
        msg = f"station_name {station_name} "
        msg += "already exists. Station names must be unique."
        msg += f"{keys}"
        assert station_name not in keys, msg
        return tuple(vertices)

    def _certify_schema(
        self,
        sampling_frequency,
        spinup,
        netcdf,
        harmonic_analysis,
    ):
        # certify sampling frequency
        msg = "Error: sampling_frequency must be None or an instance of "
        msg += f"type {timedelta}."
        assert isinstance(sampling_frequency, (type(None), timedelta)), msg
        # certify spinup
        msg = f"Error: spinup must be of type {bool}."
        assert isinstance(spinup, bool), msg
        # certify netcdf
        msg = f"Error: netcdf must be of type {bool}."
        assert isinstance(netcdf, bool), msg
        # certify harmonic_analysis
        msg = f"Error: harmonic_analysis must be of type {bool}."
        assert isinstance(harmonic_analysis, bool), msg
        return {
            'sampling_frequency': sampling_frequency,
            'spinup': spinup,
            'netcdf': netcdf,
            'harmonic_analysis': harmonic_analysis}

    @property
    def mesh(self):
        return self._mesh

    @property
    def tidal_forcing(self):
        return self._tidal_forcing

    @property
    def wind_forcing(self):
        return self._wind_forcing

    @property
    def spinup_time(self):
        return self._spinup_time

    @property
    def start_date(self):
        return self._start_date

    @property
    def end_date(self):
        return self._end_date

    @property
    def forcing_start_date(self):
        return self.start_date - self.spinup_time

    @property
    def spinup_factor(self):
        try:
            return self.__spinup_factor
        except AttributeError:
            return 1.

    @property
    def coldstart(self):
        return self.fort15('coldstart')

    @property
    def hotstart(self):
        return self.fort15('hotstart')

    @property
    def netcdf(self):
        return self._netcdf

    @property
    def container(self):
        return self._container.copy()

    @property
    def stations_output(self):
        return self.container['stations']

    @property
    def elevation_stations_output(self):
        return self.stations_output['elevation']

    @property
    def velocity_stations_output(self):
        return self.stations_output['velocity']

    @property
    def meteorological_stations_output(self):
        return self.stations_output['meteorological']

    @property
    def concentration_stations_output(self):
        return self.stations_output['concentration']

    @property
    def elevation_stations(self):
        return self.elevation_stations_output['collection']

    @property
    def velocity_stations(self):
        return self.velocity_stations_output['collection']

    @property
    def meteorological_stations(self):
        return self.meteorological_stations_output['collection']

    @property
    def concentration_stations(self):
        return self.concentration_stations_output['collection']

    @property
    def surface_outputs(self):
        return self.container['surface']

    @property
    def elevation_surface_output(self):
        return self.surface_outputs['elevation']

    @property
    def velocity_surface_output(self):
        return self.surface_outputs['velocity']

    @property
    def meteorological_surface_output(self):
        return self.surface_outputs['meteorological']

    @property
    def concentration_surface_output(self):
        return self.surface_outputs['concentration']

    @property
    def _mesh(self):
        return self.__mesh

    @property
    def _tidal_forcing(self):
        return self.__tidal_forcing

    @property
    def _wind_forcing(self):
        return self.__wind_forcing

    @property
    def _spinup_time(self):
        return self.__spinup_time

    @property
    def _start_date(self):
        return self.__start_date

    @property
    def _end_date(self):
        return self.__end_date

    @property
    def _netcdf(self):
        return self.__netcdf

    @property
    def _container(self):
        try:
            return self.__container
        except AttributeError:
            # init container
            self.__container = dict()
            # init surface outputs attributes
            schema = {
                'sampling_frequency': None,
                'spinup': False,
                'netcdf': self.netcdf,
                'harmonic_analysis': False,
                }
            self.__container['surface'] = dict()
            self.__container['surface']['elevation'] = schema.copy()
            self.__container['surface']['velocity'] = schema.copy()
            self.__container['surface']['meteorological'] = schema.copy()
            self.__container['surface']['concentration'] = schema.copy()
            # init stations output attributes
            self.__container['stations'] = dict()
            self.__container['stations']['elevation'] = schema.copy()
            self.__container['stations']['velocity'] = schema.copy()
            self.__container['stations']['meteorological'] = schema.copy()
            self.__container['stations']['concentration'] = schema.copy()
            # add a 'collections' key to hold the stations.
            self.__container['stations']['elevation'].update(
                {'collection': dict()})
            self.__container['stations']['velocity'].update(
                {'collection': dict()})
            self.__container['stations']['meteorological'].update(
                {'collection': dict()})
            self.__container['stations']['concentration'].update(
                {'collection': dict()})
            return self.__container

    @property
    def _stations_output(self):
        return self._container['stations']

    @property
    def _elevation_stations_output(self):
        return self._stations_output['elevation']

    @property
    def _velocity_stations_output(self):
        return self._stations_output['velocity']

    @property
    def _meteorological_stations_output(self):
        return self._stations_output['meteorological']

    @property
    def _concentration_stations_output(self):
        return self._stations_output['concentration']

    @property
    def _elevation_stations(self):
        return self._elevation_stations_output['collection']

    @property
    def _velocity_stations(self):
        return self._velocity_stations_output['collection']

    @property
    def _meteorological_stations(self):
        return self._meteorological_stations_output['collection']

    @property
    def _concentration_stations(self):
        return self._concentration_stations_output['collection']

    @property
    def _surface_outputs(self):
        return self._container['surface']

    @property
    def _elevation_surface_output(self):
        return self._surface_outputs['elevation']

    @property
    def _velocity_surface_output(self):
        return self._surface_outputs['velocity']

    @property
    def _meteorological_surface_output(self):
        return self._surface_outputs['meteorological']

    @property
    def _concentration_surface_output(self):
        return self._surface_outputs['concentration']

    @spinup_factor.setter
    def spinup_factor(self, spinup_factor):
        spinup_factor = float(spinup_factor)
        assert (spinup_factor >= 0 and spinup_factor <= 1.)
        self.__spinup_factor = np.abs(spinup_factor)

    @_mesh.setter
    def _mesh(self, mesh):
        assert isinstance(mesh, AdcircMesh)
        if mesh.crs is None:
            msg = 'Coordinate reference system must be set before setting up a'
            msg += ' model run.'
            raise RuntimeError(msg)
        self.__mesh = mesh

    @_tidal_forcing.setter
    def _tidal_forcing(self, tidal_forcing):
        if tidal_forcing is not None:
            assert isinstance(tidal_forcing, TidalForcing)
        self.__tidal_forcing = tidal_forcing

    @_wind_forcing.setter
    def _wind_forcing(self, wind_forcing):
        if wind_forcing is not None:
            assert isinstance(wind_forcing, WindForcing)
        self.__wind_forcing = wind_forcing

    @_spinup_time.setter
    def _spinup_time(self, spinup_time):
        if spinup_time is None:
            spinup_time = timedelta(0)
        msg = f'spinup_time must be a {timedelta} instance.'
        assert isinstance(spinup_time, timedelta), msg
        if self.tidal_forcing is not None:
            self.tidal_forcing.spinup_time = spinup_time
        self.__spinup_time = spinup_time

    @_start_date.setter
    def _start_date(self, start_date):
        if start_date is None:
            if self.wind_forcing is None:
                if self.tidal_forcing is None:
                    raise AttributeError('Error: No forcings provided.')
                else:
                    if self.tidal_forcing.start_date is None:
                        raise AttributeError('Must pass start_date argument.')
                    else:
                        start_date = self.tidal_forcing.start_date
            else:
                start_date = self.wind_forcing.start_date
        msg = 'Must pass start_date argument.'
        assert isinstance(start_date, datetime), msg
        msg = f'start_date must be a {type(datetime)} instance.'
        assert isinstance(start_date, datetime), msg
        if self.tidal_forcing is not None:
            self.tidal_forcing.start_date = start_date
        if self.wind_forcing is not None:
            self.wind_forcing.start_date = start_date
        self.__start_date = start_date

    @_end_date.setter
    def _end_date(self, end_date):
        if end_date is None:
            if self.wind_forcing is None:
                if self.tidal_forcing is None:
                    raise AttributeError('Error: No forcings provided.')
                else:
                    if self.tidal_forcing.end_date is None:
                        raise AttributeError('Must pass end_date argument.')
                    else:
                        end_date = self.tidal_forcing.end_date
            else:
                end_date = self.wind_forcing.end_date
        msg = 'Must pass end_date argument.'
        assert isinstance(end_date, datetime), msg
        msg = f'end_date must be a {type(datetime)} instance.'
        assert isinstance(end_date, datetime), msg
        if self.tidal_forcing is not None:
            self.tidal_forcing.end_date = end_date
        if self.wind_forcing is not None:
            self.wind_forcing.end_date = end_date
        self.__end_date = end_date

    @_netcdf.setter
    def _netcdf(self, netcdf):
        assert isinstance(netcdf, bool)
        self.__netcdf = netcdf
