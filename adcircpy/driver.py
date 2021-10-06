from datetime import datetime, timedelta
from functools import lru_cache
import os
import pathlib
import shutil
import subprocess
import tempfile
from typing import Any, Union

import matplotlib.pyplot as plt
import numpy as np
from psutil import cpu_count
from shapely.geometry import Point

from adcircpy.forcing import Tides  # , Winds
from adcircpy.forcing.winds.best_track import BestTrackForcing
from adcircpy.fort15 import Fort15, StationType
from adcircpy.mesh import AdcircMesh
from adcircpy.outputs.collection import OutputCollection
from adcircpy.server import SlurmConfig, SSHConfig
from adcircpy.server.driver_file import DriverFile


class AdcircRun(Fort15):
    def __init__(
        self,
        mesh: AdcircMesh,
        start_date: datetime = None,
        end_date: datetime = None,
        spinup_time: timedelta = None,
        netcdf: bool = True,
        server_config: Union[int, SSHConfig, SlurmConfig] = None,
    ):
        super().__init__(mesh)
        self._start_date = start_date
        self._end_date = end_date
        self._spinup_time = spinup_time
        self._netcdf = netcdf
        self._server_config = server_config

    def add_elevation_output_station(
        self,
        station_name: str,
        vertices: np.array
        # TODO: Is there a way to be more concise about this? (e.g. specify dimentionality?)
    ):
        self._certify_station('elevation', station_name, vertices)
        self._elevation_stations[station_name] = vertices

    def add_velocity_output_station(self, station_name: str, vertices: np.array):
        self._certify_station('velocity', station_name, vertices)
        self._velocity_stations[station_name] = vertices

    def add_meteorological_output_station(self, station_name: str, vertices: np.array):
        self._certify_station('meteorological', station_name, vertices)
        self._meteorological_stations[station_name] = vertices

    def add_concentration_output_station(self, station_name: str, vertices: np.array):
        self._certify_station('concentration', station_name, vertices)
        self._concentration_stations[station_name] = vertices

    def set_elevation_stations_output(
        self,
        sampling_rate: timedelta,
        start: datetime = None,
        end: datetime = None,
        spinup: Union[timedelta, int] = None,
        spinup_start: datetime = None,
        spinup_end: datetime = None,
        netcdf: bool = True,
        harmonic_analysis: bool = False,
    ):
        self._certify_output_request(
            sampling_rate,
            start,
            end,
            spinup,
            spinup_start,
            spinup_end,
            netcdf,
            harmonic_analysis,
        )
        self._container['stations']['elevation'].update(
            {
                'sampling_rate': sampling_rate,
                'start': start,
                'end': end,
                'spinup': spinup,
                'spinup_start': spinup_start,
                'spinup_end': spinup_end,
                'netcdf': netcdf,
                'harmonic_analysis': harmonic_analysis,
            }
        )

    def set_velocity_stations_output(
        self,
        sampling_rate: timedelta,
        start: Union[timedelta, int] = None,
        end: Union[timedelta, int] = None,
        spinup: Union[timedelta, int] = None,
        spinup_start: Union[timedelta, int] = None,
        spinup_end: Union[timedelta, int] = None,
        netcdf: bool = True,
        harmonic_analysis: bool = False,
    ):
        self._certify_output_request(
            sampling_rate,
            start,
            end,
            spinup,
            spinup_start,
            spinup_end,
            netcdf,
            harmonic_analysis,
        )
        self._container['stations']['velocity'].update(
            {
                'sampling_rate': sampling_rate,
                'start': start,
                'end': end,
                'spinup': spinup,
                'spinup_start': spinup_start,
                'spinup_end': spinup_end,
                'netcdf': netcdf,
                'harmonic_analysis': harmonic_analysis,
            }
        )

    def set_meteorological_stations_output(
        self,
        sampling_rate: timedelta,
        start: Union[timedelta, int] = None,
        end: Union[timedelta, int] = None,
        spinup: Union[timedelta, int] = None,
        spinup_start: Union[timedelta, int] = None,
        spinup_end: Union[timedelta, int] = None,
        netcdf: bool = True,
        harmonic_analysis: bool = False,
    ):
        self._certify_output_request(
            sampling_rate,
            start,
            end,
            spinup,
            spinup_start,
            spinup_end,
            netcdf,
            harmonic_analysis,
        )
        self._container['stations']['meteorological'].update(
            {
                'sampling_rate': sampling_rate,
                'start': start,
                'end': end,
                'spinup': spinup,
                'spinup_start': spinup_start,
                'spinup_end': spinup_end,
                'netcdf': netcdf,
                'harmonic_analysis': harmonic_analysis,
            }
        )

    def set_concentration_stations_output(
        self,
        sampling_rate: timedelta,
        start: Union[timedelta, int] = None,
        end: Union[timedelta, int] = None,
        spinup: Union[timedelta, int] = None,
        spinup_start: Union[timedelta, int] = None,
        spinup_end: Union[timedelta, int] = None,
        netcdf: bool = True,
        harmonic_analysis: bool = False,
    ):
        self._certify_output_request(
            sampling_rate,
            start,
            end,
            spinup,
            spinup_start,
            spinup_end,
            netcdf,
            harmonic_analysis,
        )
        self._container['stations']['concentration'].update(
            {
                'sampling_rate': sampling_rate,
                'start': start,
                'end': end,
                'spinup': spinup,
                'spinup_start': spinup_start,
                'spinup_end': spinup_end,
                'netcdf': netcdf,
                'harmonic_analysis': harmonic_analysis,
            }
        )

    def set_elevation_surface_output(
        self,
        sampling_rate: timedelta,
        start: Union[timedelta, int] = None,
        end: Union[timedelta, int] = None,
        spinup: Union[timedelta, int] = None,
        spinup_start: Union[timedelta, int] = None,
        spinup_end: Union[timedelta, int] = None,
        netcdf: bool = True,
        harmonic_analysis: bool = False,
    ):
        self._certify_output_request(
            sampling_rate,
            start,
            end,
            spinup,
            spinup_start,
            spinup_end,
            netcdf,
            harmonic_analysis,
        )
        self._container['surface']['elevation'].update(
            {
                'sampling_rate': sampling_rate,
                'start': start,
                'end': end,
                'spinup': spinup,
                'spinup_start': spinup_start,
                'spinup_end': spinup_end,
                'netcdf': netcdf,
                'harmonic_analysis': harmonic_analysis,
            }
        )

    def set_velocity_surface_output(
        self,
        sampling_rate: timedelta,
        start: Union[timedelta, int] = None,
        end: Union[timedelta, int] = None,
        spinup: Union[timedelta, int] = None,
        spinup_start: Union[timedelta, int] = None,
        spinup_end: Union[timedelta, int] = None,
        netcdf: bool = True,
        harmonic_analysis: bool = False,
    ):
        self._certify_output_request(
            sampling_rate,
            start,
            end,
            spinup,
            spinup_start,
            spinup_end,
            netcdf,
            harmonic_analysis,
        )
        self._container['surface']['velocity'].update(
            {
                'sampling_rate': sampling_rate,
                'start': start,
                'end': end,
                'spinup': spinup,
                'spinup_start': spinup_start,
                'spinup_end': spinup_end,
                'netcdf': netcdf,
                'harmonic_analysis': harmonic_analysis,
            }
        )

    def set_meteorological_surface_output(
        self,
        sampling_rate: timedelta,
        start: Union[timedelta, int] = None,
        end: Union[timedelta, int] = None,
        spinup: Union[timedelta, int] = None,
        spinup_start: Union[timedelta, int] = None,
        spinup_end: Union[timedelta, int] = None,
        netcdf: bool = True,
        harmonic_analysis: bool = False,
    ):
        self._certify_output_request(
            sampling_rate,
            start,
            end,
            spinup,
            spinup_start,
            spinup_end,
            netcdf,
            harmonic_analysis,
        )
        self._container['surface']['meteorological'].update(
            {
                'sampling_rate': sampling_rate,
                'start': start,
                'end': end,
                'spinup': spinup,
                'spinup_start': spinup_start,
                'spinup_end': spinup_end,
                'netcdf': netcdf,
                'harmonic_analysis': harmonic_analysis,
            }
        )

    def set_concentration_surface_output(
        self,
        sampling_rate: timedelta,
        start: Union[timedelta, int] = None,
        end: Union[timedelta, int] = None,
        spinup: Union[timedelta, int] = None,
        spinup_start: Union[timedelta, int] = None,
        spinup_end: Union[timedelta, int] = None,
        netcdf: bool = True,
        harmonic_analysis: bool = False,
    ):
        self._certify_output_request(
            sampling_rate,
            start,
            end,
            spinup,
            spinup_start,
            spinup_end,
            netcdf,
            harmonic_analysis,
        )
        self._container['surface']['concentration'].update(
            {
                'sampling_rate': sampling_rate,
                'start': start,
                'end': end,
                'spinup': spinup,
                'spinup_start': spinup_start,
                'spinup_end': spinup_end,
                'netcdf': netcdf,
                'harmonic_analysis': harmonic_analysis,
            }
        )

    def remove_elevation_output_station(self, station_name):
        self._elevation_stations.pop(station_name)

    def remove_velocity_output_station(self, station_name):
        self._velocity_stations.pop(station_name)

    def remove_meteorological_output_station(self, station_name):
        self._meteorological_stations.pop(station_name)

    def remove_concentration_output_station(self, station_name):
        self._concentration_stations.pop(station_name)

    def write(
        self,
        output_directory: str,
        overwrite: bool = False,
        fort14: str = 'fort.14',
        fort13: str = 'fort.13',
        fort22: str = 'fort.22',
        fort15: str = 'fort.15',
        coldstart: str = 'fort.15.coldstart',
        hotstart: str = 'fort.15.hotstart',
        driver: str = 'driver.sh',
        nproc: int = None,
    ):
        output_directory = pathlib.Path(output_directory)
        output_directory.mkdir(parents=True, exist_ok=True)

        # write fort.14
        if fort14:
            path = output_directory / fort14
            self.mesh.write(path, overwrite)

        # write fort.13 (optional)
        if fort13:
            if len(self.mesh.get_nodal_attribute_names()) > 0:
                self.mesh.nodal_attributes.write(output_directory / fort13, overwrite)

        # write fort.15 -> this part has two different cases that depend
        # on the input configuration given by the user.

        # CASE 1:
        # When wind forcing is present, the run has to (almost) necessarily
        # be executed in two phases: coldstart and hotstart.
        # In a nutshell, a single phase run (coldstart only) that includes
        # meteorological inputs is currently not supported.
        # For this reason, no spiunp time given means a single phase run
        # which at the current stage implies tides only.
        # In this case we set IHOT=0 but call the hotstart writer.
        if self.spinup_time == timedelta(seconds=0):
            # easiest way is to override IHOT to 0
            # IHOT depends on _runtype which is not set on this case.
            self._IHOT = 0
            # and call the hotstart writer,
            super().write('hotstart', output_directory / fort15, overwrite)
            if self.wind_forcing is not None:
                if fort22:
                    self.wind_forcing.write(output_directory / fort22, overwrite)

        # CASE 2:
        # This is a run that the user specified some ramping time.
        # It may or may not include wind files.
        # In other words, one can do a tidal only run as coldstart/hotstart
        # as well as a tidal-only spinup + meteorological hotstart.
        else:
            if self.wind_forcing is not None:
                if fort22:
                    self.wind_forcing.write(output_directory / fort22, overwrite)
            if coldstart:
                super().write('coldstart', output_directory / coldstart, overwrite)
            if hotstart:
                super().write('hotstart', output_directory / hotstart, overwrite)

        if driver is not None:
            if isinstance(self._server_config, SlurmConfig):
                driver = self._server_config._filename
            script = DriverFile(self, nproc)
            script.write(output_directory / driver, overwrite)

    def import_stations(
        self,
        fort15: os.PathLike,
        station_types: [StationType] = None,
        only_within: bool = False,
    ):
        if station_types is None:
            station_types = [StationType.ELEVATION, StationType.VELOCITY]
            if self.IM == 10:
                station_types.append(StationType.CONCENTRATION)
            if self.NWS > 0:
                station_types.append(StationType.METEOROLOGICAL)

        envelope = self.mesh.hull.rings.multipolygon if only_within else None
        stations = Fort15.parse_stations(path=fort15, station_types=station_types)
        for station_type, station_vertices in stations.items():
            for name, vertex in station_vertices.items():
                if not only_within or Point(vertex).within(envelope):
                    if station_type == StationType.ELEVATION:
                        self.add_elevation_output_station(name, vertex)
                    elif station_type == StationType.VELOCITY:
                        self.add_velocity_output_station(name, vertex)
                    elif station_type == StationType.CONCENTRATION:
                        self.add_concentration_output_station(name, vertex)
                    elif station_type == StationType.METEOROLOGICAL:
                        self.add_meteorological_output_station(name, vertex)

    def run(
        self,
        outdir=None,
        nproc=-1,
        overwrite=False,
        coldstart=True,
        hotstart=True,
        server_config=None,
    ):

        if outdir is None:
            self._outdir_tmpdir = tempfile.TemporaryDirectory()
            outdir = pathlib.Path(self._outdir_tmpdir.name)
        else:
            outdir = pathlib.Path(outdir)
        #     if outdir.exists() and not overwrite:
        #         msg = f"{outdir} exists and overwrite is not enabled."
        #         raise IOError(msg)

        # if not outdir.exists():
        #     outdir.mkdir(parents=True)

        # local adcirc run

        if server_config is None:
            self._run_local(
                nproc=nproc,
                outdir=outdir,
                overwrite=overwrite,
                coldstart=coldstart,
                hotstart=hotstart,
            )

            # server adcirc run
        else:
            server_config.run(
                driver=self,
                outdir=outdir,
                overwrite=overwrite,
                coldstart=coldstart,
                hotstart=hotstart,
            )

        self._load_outdir(outdir)

        return self._output_collection

    @property
    def mesh(self):
        return self._mesh

    @property
    def tidal_forcing(self):
        if not hasattr(self, '_tidal_forcing'):
            self._tidal_forcing = self.mesh.forcings.tides
            if isinstance(self._tidal_forcing, Tides):
                self._tidal_forcing.start_date = self.start_date
                self._tidal_forcing.end_date = self.end_date
                self._tidal_forcing.spinup_time = self.spinup_time
        return self._tidal_forcing

    @property
    def wind_forcing(self):
        return self.mesh.forcings.wind

    @property
    def wave_forcing(self):
        return self.mesh.forcings.wave

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
            return 1.0

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
    def output_collection(self):
        return self._output_collection

    def _load_outdir(self, outdir):

        # gather outputs
        maxele = pathlib.Path(outdir / 'hotstart/maxele.63.nc')

        self._output_collection = OutputCollection(
            maxele=maxele if maxele.is_file() else None, crs=self.mesh.crs
        )

        if len(self._output_collection) == 0:
            raise Exception('No outputs found.')

        # for output in self._output_collection:
        #     output.make_plot(show=True)

        return self.output_collection

    def _run_padcirc(self, rundir, nproc=-1):
        nproc = self._get_nproc(nproc)
        cmd = list()
        cmd.append('mpiexec')
        cmd.append('-n')
        cmd.append(str(nproc))
        cmd.append('padcirc')
        err = self._launch_command(cmd, rundir)

        msg = '** ERROR: Elevation.gt.ErrorElev, ADCIRC stopping. **'
        if msg in "".join(err):
            print(msg)
            self._handle_blowup(err)

        # filter IEEE_UNDERFLOW_FLAG IEEE_DENORMAL
        msg = 'Note: The following floating-point exceptions are signalling:'
        err = [line for line in err if msg not in line]
        if len(err) > 0:
            if msg not in "".join(err):
                msg = '\n'
                msg += "".join(err)
                raise Exception(msg)
            else:
                raise Exception(msg)

        return err

    def _run_local(self, nproc, outdir, overwrite, coldstart, hotstart):
        self.write(outdir, overwrite, driver=None)
        if self.spinup_time > timedelta(seconds=0):
            if coldstart:
                self._run_coldstart(nproc, outdir)
            if hotstart:
                self._run_hotstart(nproc, outdir)
        else:
            # not separated into coldstart/hotstart
            self._run_single_phase(nproc, outdir)

    def _run_coldstart(self, nproc, wdir):
        self._stage_files('coldstart', nproc, wdir)
        self._run_adcprep('coldstart', nproc, wdir)
        self._run_padcirc(wdir / 'coldstart', nproc)

    def _run_hotstart(self, nproc, wdir):
        self._stage_files('hotstart', nproc, wdir)
        self._run_adcprep('hotstart', nproc, wdir)
        if isinstance(self.mesh.forcings.wind, BestTrackForcing):
            self._run_aswip(wdir / 'hotstart')
        if self.wave_forcing is not None:
            if self.waves_focing.model.lower() == 'swan':
                self._run_padcswan(wdir / 'hotstart', nproc)
            else:
                msg = 'Unknown wave coupling type.'
                raise NotImplementedError(msg)
        else:
            self._run_padcirc(wdir / 'hotstart', nproc)

    def _run_aswip(self, wdir):
        subprocess.check_call(['aswip'], cwd=wdir)
        (wdir / f'NWS_{self.mesh.forcings.wind.NWS}_fort.22').replace(wdir / 'fort.22')

    def _run_single_phase(self, nproc, wdir):
        # "single phase" means not separated into coldstart/hotstart
        self._run_adcprep('./', nproc, wdir)
        self._run_padcirc(wdir, nproc)

    def _stage_files(self, runtype, nproc, wdir):

        # create run directory
        cwdir = wdir / runtype
        if cwdir.exists():
            shutil.rmtree(cwdir.absolute())
        cwdir.mkdir()

        # symlink files
        os.symlink(wdir / 'fort.14', cwdir / 'fort.14')
        os.symlink(wdir / 'fort.13', cwdir / 'fort.13')
        os.symlink(wdir / f'fort.15.{runtype}', cwdir / 'fort.15')
        if isinstance(self.mesh.forcings.wind, BestTrackForcing):
            os.symlink(wdir / 'fort.22', cwdir / 'fort.22')
        if runtype == 'hotstart':
            os.symlink(wdir / 'coldstart/fort.67.nc', cwdir / 'fort.67.nc')

    def _run_adcprep(self, runtype, nproc, wdir):

        # run directory
        cwdir = wdir / runtype
        nproc = self._get_nproc(nproc)

        if nproc > 1:
            cmd = list()
            cmd.append('adcprep')
            cmd.append('--np')
            cmd.append(f'{nproc:d}')
            self._launch_command(cmd + ['--partmesh'], cwdir)
            self._launch_command(cmd + ['--prepall'], cwdir)

        else:
            raise NotImplementedError('run serial')

    def _certify_station(self, physical_var, station_name, vertices):
        keys = list(self._container['stations'][physical_var].keys())
        vertices = np.asarray(vertices)
        msg = 'vertices argument must be a two-tuple of floats.'
        assert vertices.shape == (2,), msg
        msg = f'station_name {station_name} '
        msg += 'already exists. Station names must be unique.'
        msg += f'{keys}'
        assert station_name not in keys, msg
        return tuple(vertices)

    def _handle_blowup(self, err):
        data = self._get_blowup_data(err)
        # TODOL: There's an ambiguety here, we don't know if ADCIRC is
        # reporting the FORTRAN index or the node hash. We assume the former
        # case here.
        idx = np.asarray(data['maxele_node']) - 1
        ax = self.mesh.make_plot()
        self.mesh.triplot(axes=ax)
        ax.scatter(self.mesh.x[idx], self.mesh.y[idx], s=80, facecolors='none', edgecolors='r')
        ax.axis('scaled')
        plt.show()

    def _certify_output_request(
        self,
        sampling_rate: timedelta,
        start: datetime,
        end: datetime,
        spinup: Union[timedelta, int],
        spinup_start: datetime,
        spinup_end: datetime,
        netcdf: bool,
        harmonic_analysis: bool,
    ):
        self._validate_argument(sampling_rate, timedelta, 'sampling_rate')
        self._validate_argument(start, datetime, 'start')
        self._validate_argument(end, datetime, 'end')
        self._validate_argument(spinup, [timedelta, int], 'spinup')
        self._validate_argument(spinup_start, datetime, 'spinup_start')
        self._validate_argument(spinup_end, datetime, 'spinup_end')
        self._validate_argument(netcdf, bool, 'netcdf', include_none=False)
        self._validate_argument(
            harmonic_analysis, bool, 'harmonic_analysis', include_none=False
        )

    def _write_bash_driver(self, destination):
        source = pathlib.Path(__file__).parent / 'padcirc_driver.sh'
        shutil.copyfile(source, destination)

    @staticmethod
    def _validate_argument(
        value: Any, types: [type], name: str = None, include_none: bool = True
    ):
        if isinstance(types, type):
            types = [types]
        if include_none:
            types.append(type(None))
        types = tuple(types)

        name = f'"{name}"' if name is not None else 'value'

        if not isinstance(value, types):
            raise TypeError(f'{name} is not of type(s) {types}')

    @staticmethod
    def _launch_command(cmd, rundir):
        p = subprocess.Popen(
            cmd,
            universal_newlines=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=rundir.absolute(),
        )
        for line in p.stdout:
            if 'MPI terminated with Status =' in line:
                break
            print(line, end="")
        p.wait()
        lines = p.stderr.readlines()
        # p.close()
        return lines

    @staticmethod
    def _get_nproc(nproc):
        return cpu_count(logical=False) if nproc == -1 else nproc

    @staticmethod
    def _get_blowup_data(err):
        s = "".join(err).split('** WARNING: Elevation.gt.WarnElev **')
        s = [_ for _ in s if 'TIME' in _]
        blowup = {
            'timestep': list(),
            'time': list(),
            'maxele': list(),
            'maxele_node': list(),
            'maxvel': list(),
            'maxvel_node': list(),
        }
        for output in s:
            blowup['timestep'].append(int(output.split('TIME STEP =')[1].split()[0]))
            blowup['time'].append(float(output.split('TIME =')[1].split()[0]))
            blowup['maxele'].append(float(output.split('ELMAX =')[1].split()[0]))
            blowup['maxvel'].append(float(output.split('SPEEDMAX =')[1].split()[0]))
            nodes = output.split('AT NODE')
            blowup['maxele_node'].append(int(nodes[1].split()[0]))
            blowup['maxvel_node'].append(int(nodes[2].split()[0]))
        return blowup

    @property
    def _mesh(self):
        return self.__mesh

    @_mesh.setter
    def _mesh(self, mesh):
        assert isinstance(mesh, AdcircMesh)
        if mesh.crs is None:
            msg = 'Coordinate reference system must be set before setting up a'
            msg += ' model run.'
            raise RuntimeError(msg)
        self.__mesh = mesh

    @property
    def _spinup_time(self):
        return self.__spinup_time

    @_spinup_time.setter
    def _spinup_time(self, spinup_time):
        if spinup_time is None:
            spinup_time = timedelta(0)
        assert isinstance(spinup_time, timedelta)
        self.__spinup_time = spinup_time

    @property
    def _start_date(self):
        return self.__start_date

    @_start_date.setter
    def _start_date(self, start_date):
        if start_date is None:
            if isinstance(self.mesh.forcings.wind, BestTrackForcing):
                start_date = self.mesh.forcings.wind.start_date
                self.mesh.forcings.tides.start_date = self.mesh.forcings.wind.start_date
        else:
            if isinstance(self.mesh.forcings.tides, Tides):
                self.mesh.forcings.tides.start_date = start_date
            if isinstance(self.mesh.forcings.wind, BestTrackForcing):
                self.mesh.forcings.wind.start_date = start_date
        assert isinstance(start_date, datetime)
        self.__start_date = start_date

    @property
    def _end_date(self):
        return self.__end_date

    @_end_date.setter
    def _end_date(self, end_date):
        if isinstance(end_date, timedelta):
            end_date = self._start_date + end_date
        if end_date is None:
            if isinstance(self.wind_forcing, BestTrackForcing):
                end_date = self.wind_forcing.end_date
                self.mesh.forcings.tides.end_date = self.mesh.forcings.wind.end_date
        else:
            if isinstance(self.mesh.forcings.tides, Tides):
                self.mesh.forcings.tides.end_date = end_date
            if isinstance(self.mesh.forcings.wind, BestTrackForcing):
                self.mesh.forcings.wind.end_date = end_date

        assert isinstance(end_date, datetime)
        assert end_date > self.start_date
        self.__end_date = end_date

    @property
    def _netcdf(self):
        return self.__netcdf

    @_netcdf.setter
    def _netcdf(self, netcdf):
        assert isinstance(netcdf, bool)
        self.__netcdf = netcdf

    @property
    def _server_config(self):
        return self.__server_config

    @_server_config.setter
    def _server_config(self, server_config):
        if server_config is None:
            server_config = self._get_nproc(-1)
        msg = 'server_config must be int, SSHConfig or SlurmConfig'
        assert isinstance(server_config, (int, SSHConfig, SlurmConfig)), msg
        self.__server_config = server_config

    @property
    @lru_cache(maxsize=None)
    def _container(self):
        # init container
        container = dict()
        # init surface outputs attributes
        schema = {
            'sampling_rate': None,
            'start': None,
            'end': None,
            'spinup': None,
            'spinup_start': None,
            'spinup_end': None,
            'netcdf': self.netcdf,
            'harmonic_analysis': False,
        }

        for otype in ['surface', 'stations']:
            container[otype] = dict()
            for ovar in ['elevation', 'velocity', 'meteorological', 'concentration']:
                container[otype][ovar] = schema.copy()
                if otype == 'stations':
                    container[otype][ovar].update({'collection': dict()})
        return container

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

    @property
    def _output_collection(self):
        return self.__output_collection

    @_output_collection.setter
    def _output_collection(self, output_collection):
        assert isinstance(output_collection, OutputCollection)
        self.__output_collection = output_collection
