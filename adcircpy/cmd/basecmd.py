# import abc
from datetime import timedelta
from functools import lru_cache
import pathlib

from adcircpy import AdcircMesh, AdcircRun, Tides, server


class AdcircCommand:

    def __init__(self, args):
        self._args = args

    def run(self):

        # write and exit if generate only
        if self._args.generate_only:
            self.driver.write(
                self._args.output_directory,
                overwrite=self._args.overwrite
            )
            return

        outputs = self.driver.run(
            outdir=self.output_directory,
            nproc=self._args.nproc,
            overwrite=self._args.overwrite,
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
            server_config=self.server_config
        )
        self._enable_outputs(driver)
        if self._args.timestep:
            driver.timestep = self._args.timestep
        driver.gwce_solution_scheme = self._args.gwce_solution_scheme
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
            self.__spinup_time = timedelta(days=self._args.spinup_days)
            return self.__spinup_time

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
        if self._args.output_directory is not None:
            return pathlib.Path(self._args.output_directory).absolute()

    @property
    def constituents(self):
        try:
            return self.__constituents
        except AttributeError:
            # might be better to get these from Tides()
            _major = ('Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2')
            _all = (*_major, 'Mm', 'Mf', 'M4', 'MN4', 'MS4', '2N2', 'S1')
            if ('all' in self._args.constituents
                    and len(self._args.constituents) > 1):
                msg = 'When using all, must only pass one'
                raise IOError(msg)

            elif ('major' in self._args.constituents
                  and len(self._args.constituents) > 1):
                msg = 'When using major, must only pass one'
                raise IOError(msg)
            if 'all' in self._args.constituents:
                constituents = _all
            elif 'major' in self._args.constituents:
                constituents = _major
            else:
                constituents = self._args.constituents
            self.__constituents = constituents
            return self.__constituents

    @property
    @lru_cache(maxsize=None)
    def server_config(self):
        if self._args.hostname:
            if (not self._args.use_slurm or
                    not self._args.use_torque or
                    not self._args.use_pbs):
                return server.ServerConfig(
                    hostname=self._args.hostname,
                    nprocs=self._args.nproc,
                    wdir=self._args.wdir,
                    binaries_prefix=self._args.binaries_prefix,
                    source_script=self._args.source_script,
                    additional_mpi_options=self._args.additional_mpi_options,
                )

        if self._args.use_slurm:
            kwargs = {
                "account": self._args.account,
                "ntasks": self._args.slurm_ntasks,
                "partition": self._args.partition,
                "walltime": timedelta(hours=self._args.walltime),
                "mail_type": self._args.mail_type,
                "mail_user": self._args.mail_user,
                "log_filename": self._args.log_filename,
                "modules": self._args.modules,
                "path_prefix": self._args.path_prefix,
                "extra_commands": self._args.extra_commands,
                "launcher": self._args.slurm_launcher,
                "nodes": self._args.slurm_nodes
            }
            if self._args.slurm_filename is not None:
                kwargs.update({"filename": self._args.slurm_ntasks})
            if self._args.slurm_rundir is not None:
                kwargs.update({"run_directory": self._args.slurm_rundir})
            if self._args.run_name is not None:
                kwargs.update({"run_name": self._args.run_name})

            return server.SlurmConfig(**kwargs)

        elif self._args.use_torque or self._args.use_pbs:
            raise NotImplementedError

    def _enable_outputs(self, driver):
        self._enable_output(driver, 'elevation', 'surface')
        self._enable_output(driver, 'velocity', 'surface')
        self._enable_output(driver, 'meteorological', 'surface')
        self._enable_output(driver, 'concentration', 'surface')
        self._init_output_stations(driver)

    def _enable_output(self, driver, name, _type):
        fs = getattr(self._args, f"{name}_{_type}_sampling_rate")
        if fs is not None:
            fs = timedelta(minutes=fs)
        fss = getattr(self._args, f"{name}_{_type}_sampling_rate_spinup")
        if fss is not None:
            fss = timedelta(minutes=fss)
        ha = getattr(self._args, f"{name}_{_type}_harmonic_analysis")
        # has = getattr(self._args, f"{name}_{_type}_harmonic_analysis_spinup")
        getattr(driver, f"set_{name}_{_type}_output")(
            sampling_rate=fs,
            harmonic_analysis=ha,
            spinup=fss,
            netcdf=self._args.netcdf,
        )

    def _init_output_stations(self, driver):
        if self._args.stations_file is not None:
            driver.import_stations(
                pathlib.Path(self._args.stations_file).resolve())
            self._enable_output(driver, 'elevation', 'stations')
            self._enable_output(driver, 'velocity', 'stations')
            self._enable_output(driver, 'meteorological', 'stations')
            self._enable_output(driver, 'concentration', 'stations')

    @property
    @lru_cache(maxsize=None)
    def _mesh(self):
        mesh = AdcircMesh.open(
            self._args.mesh,
            self._args.crs
        )

        if self._args.generate_boundaries:
            mesh.generate_boundaries(
                threshold=self._args.boundaries_threshold,
                land_ibtype=self._args.land_ibtype,
                interior_ibtype=self._args.island_ibtype,
            )

        # set nodal attributes
        if self._args.fort13 is not None:
            mesh.import_nodal_attributes(
                pathlib.Path(self._args.fort13).resolve()
            )

        if 'all' in self._args.coldstart_attributes:
            for attr in mesh.get_nodal_attribute_names():
                mesh.set_nodal_attribute_coldstart_state(attr, True)
        else:
            for attr in self._args.coldstart_attributes:
                mesh.set_nodal_attribute_coldstart_state(attr, True)

        if 'all' in self._args.hotstart_attributes:
            for attr in mesh.get_nodal_attribute_names():
                mesh.set_nodal_attribute_hotstart_state(attr, True)
        else:
            for attr in self._args.hotstart_attributes:
                mesh.set_nodal_attribute_hotstart_state(attr, True)

        if self._args.generate_tau0:
            mesh.generate_tau0()

        if self.tidal_forcing is not None:
            mesh.add_forcing(self.tidal_forcing)

        if self.wave_forcing is not None:
            mesh.add_forcing(self.wave_forcing)

        if self.wind_forcing is not None:
            mesh.add_forcing(self.wind_forcing)
        return mesh
