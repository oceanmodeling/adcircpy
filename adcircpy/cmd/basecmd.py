# import abc
from datetime import timedelta
from functools import lru_cache
import pathlib

from adcircpy import AdcircMesh, AdcircRun, Tides  # BestTrackForcing, ServerConfig, SlurmConfig,


class AdcircCommand:

    def __init__(self, args):
        self._args = args

    def run(self):

        # write and exit if generate only
        if self.args.generate_only:
            self.driver.write(
                self.args.output_directory,
                overwrite=self.args.overwrite
                )
            return

        outputs = self.driver.run(
            outdir=self.output_directory,
            nproc=self.args.nproc,
            overwrite=self.args.overwrite,
            server_config=self.server_config,
        )
        self._output_collection = outputs

    @property
    @lru_cache(maxsize=None)
    def driver(self):
        driver = AdcircRun(
            self.mesh,
            self.start_date,
            self.end_date,
            self.spinup_time,
        )
        self._enable_outputs(driver)
        if self.args.timestep:
            driver.timestep = self.args.timestep
        driver.gwce_solution_scheme = self.args.gwce_solution_scheme
        return driver

    @property
    def start_date(self):
        return self._start_date

    @property
    def end_date(self):
        return self._end_date

    @property
    def spinup_time(self):
        try:
            return self.__spinup_time
        except AttributeError:
            self.__spinup_time = timedelta(days=self.args.spinup_days)
            return self.__spinup_time

    @property
    def args(self):
        return self._args

    @property
    def mesh(self):
        return self._mesh

    @property
    @lru_cache(maxsize=None)
    def tidal_forcing(self):
        tidal_forcing = Tides()
        for constituent in self.constituents:
            tidal_forcing.use_constituent(constituent)
        return tidal_forcing

    @property
    def wind_forcing(self):
        try:
            return self._wind_forcing
        except AttributeError:
            return None

    @property
    def wave_forcing(self):
        try:
            return self._wave_forcing
        except AttributeError:
            return None

    @property
    def output_directory(self):
        if self.args.output_directory is not None:
            return pathlib.Path(self.args.output_directory).absolute()

    @property
    def constituents(self):
        try:
            return self.__constituents
        except AttributeError:
            # might be better to get these from Tides()
            _major = ('Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2')
            _all = (*_major, 'Mm', 'Mf', 'M4', 'MN4', 'MS4', '2N2', 'S1')
            if ('all' in self.args.constituents
                    and len(self.args.constituents) > 1):
                msg = 'When using all, must only pass one'
                raise IOError(msg)

            elif ('major' in self.args.constituents
                  and len(self.args.constituents) > 1):
                msg = 'When using major, must only pass one'
                raise IOError(msg)
            if 'all' in self.args.constituents:
                constituents = _all
            elif 'major' in self.args.constituents:
                constituents = _major
            else:
                constituents = self.args.constituents
            self.__constituents = constituents
            return self.__constituents

    @property
    def server_config(self):
        if self.args.hostname:
            if (not self.args.use_slurm or
                    not self.args.use_torque or
                    not self.args.use_pbs):
                server_config = ServerConfig(
                    hostname=self.args.hostname,
                    nprocs=self.args.nproc,
                    wdir=self.args.wdir,
                    binaries_prefix=self.args.binaries_prefix,
                    source_script=self.args.source_script,
                    additional_mpi_options=self.args.additional_mpi_options,
                )

            elif self.args.use_slurm:
                raise NotImplementedError
                server_config = SlurmConfig(
                )

            elif self.args.use_torque or self.args.use_pbs:
                raise NotImplementedError

        else:
            server_config = None

        self.__server_config = server_config
        return self.__server_config

    def _enable_outputs(self, driver):
        self._enable_output(driver, 'elevation', 'surface')
        self._enable_output(driver, 'velocity', 'surface')
        self._enable_output(driver, 'meteorological', 'surface')
        self._enable_output(driver, 'concentration', 'surface')
        self._init_output_stations(driver)

    def _enable_output(self, driver, name, _type):
        fs = getattr(self.args, f"{name}_{_type}_sampling_rate")
        if fs is not None:
            fs = timedelta(minutes=fs)
        fss = getattr(self.args, f"{name}_{_type}_sampling_rate_spinup")
        if fss is not None:
            fss = timedelta(minutes=fss)
        ha = getattr(self.args, f"{name}_{_type}_harmonic_analysis")
        # has = getattr(self.args, f"{name}_{_type}_harmonic_analysis_spinup")
        getattr(driver, f"set_{name}_{_type}_output")(
            sampling_rate=fs,
            harmonic_analysis=ha,
            spinup=fss,
            netcdf=self.args.netcdf,
        )

    def _init_output_stations(self, driver):
        if self.args.stations_file is not None:
            driver.import_stations(
                pathlib.Path(self.args.stations_file).resolve())
            self._enable_output(driver, 'elevation', 'stations')
            self._enable_output(driver, 'velocity', 'stations')
            self._enable_output(driver, 'meteorological', 'stations')
            self._enable_output(driver, 'concentration', 'stations')

    @property
    def _args(self):
        return self.__args

    @property
    @lru_cache(maxsize=None)
    def _mesh(self):
        mesh = AdcircMesh.open(
            self.args.mesh,
            self.args.crs
        )

        if self.args.generate_boundaries:
            mesh.generate_boundaries(
                threshold=self.args.boundaries_threshold,
                land_ibtype=self.args.land_ibtype,
                interior_ibtype=self.args.island_ibtype,
                )

        # set nodal attributes
        if self.args.fort13 is not None:
            mesh.import_nodal_attributes(
                pathlib.Path(self.args.fort13).resolve()
            )

        if 'all' in self.args.coldstart_attributes:
            for attr in mesh.get_nodal_attribute_names():
                mesh.set_nodal_attribute_coldstart_state(attr, True)
        else:
            for attr in self.args.coldstart_attributes:
                mesh.set_nodal_attribute_coldstart_state(attr, True)

        if 'all' in self.args.hotstart_attributes:
            for attr in mesh.get_nodal_attribute_names():
                mesh.set_nodal_attribute_hotstart_state(attr, True)
        else:
            for attr in self.args.hotstart_attributes:
                mesh.set_nodal_attribute_hotstart_state(attr, True)

        if self.args.generate_tau0:
            mesh.generate_tau0()

        if self.tidal_forcing is not None:
            mesh.add_forcing(self.tidal_forcing)

        if self.wave_forcing is not None:
            mesh.add_forcing(self.wave_forcing)

        if self.wind_forcing is not None:
            mesh.add_forcing(self.wind_forcing)
        return mesh

    @_args.setter
    def _args(self, args):
        self.__args = args
