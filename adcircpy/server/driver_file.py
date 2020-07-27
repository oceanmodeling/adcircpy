import pathlib


from adcircpy.server.slurm import SlurmConfig


class DriverFile:

    def __init__(self, driver):
        self._driver = driver

    def write(self, path, overwrite=False):
        with open(pathlib.Path(path), 'w') as f:
            f.write(self._script)
        # os.chmod(path, 744)

    @property
    def _script(self):
        f = f"{self._shebang}\n"
        if not isinstance(self._server_config, int):
            f += self._server_config._prefix
        else:
            f += '\nset -e\n'

        if self._executable.startswith('p') and isinstance(self._server_config, int):
            if self._server_config > 1:
                f += f"\nNPROCS={self._nprocs}\n"

        if self._driver.spinup_time.total_seconds() == 0:
            f += self._single_phase_run

        else:
            f += self._dual_phase_run

        f += self._clean_directory

        f += "\nmain\n"

        return f

    @property
    def _shebang(self):
        f = '#!/bin/bash'
        if not isinstance(self._server_config, int):
            f += ' --login'
        return f

    @property
    def _single_phase_run(self):
        f = r"""
main()
{
  rm -rf work
  mkdir work
  cd work
  ln -sf ../fort.14
  ln -sf ../fort.13
  ln -sf ../fort.15 ./fort.15
"""

        if self._executable.startswith('p'):
            f += f"  adcprep --np {self._nprocs} --partmesh\n"
            f += f"  adcprep --np {self._nprocs} --prepall\n"
            f += f"  {self._mpi} {self._executable}"
        else:
            f += f"  {self._executable}"
        f += """
  rm -rf fort.68.nc
  clean_directory
  cd ..
}
"""
        return f

    @property
    def _dual_phase_run(self):
        f = ''
        f += self._bash_main_dual_phase
        f += self._run_coldstart_phase
        f += self._run_hotstart_phase
        return f

    @property
    def _bash_main_dual_phase(self):
        f = """
main() {
  SECONDS=0
  run_coldstart_phase
  """
        f += f'if grep -Rq "ERROR: Elevation.gt.ErrorElev, ADCIRC stopping." {self._logfile}; then'
        f += """
    duration=$SECONDS
    echo "ERROR: Elevation.gt.ErrorElev, ADCIRC stopping."
    echo "Wallclock time: $(($duration / 60)) minutes and $(($duration % 60)) seconds."
    exit -1
  else
    run_hotstart_phase
    duration=$SECONDS
    """
        f += f'if grep -Rq "ERROR: Elevation.gt.ErrorElev, ADCIRC stopping." {self._logfile}; then'
        f += """
      echo "ERROR: Elevation.gt.ErrorElev, ADCIRC stopping."
      echo "Wallclock time: $(($duration / 60)) minutes and $(($duration % 60)) seconds."
      exit -1
    fi
  fi
  echo "Wallclock time: $(($duration / 60)) minutes and $(($duration % 60)) seconds."
}
"""
        return f

    @property
    def _run_coldstart_phase(self):
        f = """
run_coldstart_phase() {
  rm -rf coldstart
  mkdir coldstart
  cd coldstart
  ln -sf ../fort.14
  ln -sf ../fort.13
  ln -sf ../fort.15.coldstart ./fort.15
"""

        if self._executable.startswith('p'):
            if isinstance(self._server_config, SlurmConfig):
                f += "  adcprep --np $SLURM_NTASKS --partmesh\n"
                f += "  adcprep --np $SLURM_NTASKS --prepall\n"
            else:
                f += f"  adcprep --np {self._nprocs} --partmesh\n"
                f += f"  adcprep --np {self._nprocs} --prepall\n"
            f += f"  {self._mpi} {self._executable} "
        else:
            f += f"  {self._executable} "
        if not isinstance(self._server_config, SlurmConfig):
            f += f"2>&1 | tee ../{self._logfile}"
        f += """
  rm -rf fort.68.nc
  clean_directory
  cd ..
}
"""
        return f

    @property
    def _run_hotstart_phase(self):
        f = """
run_hotstart_phase() {
  rm -rf hotstart
  mkdir hotstart
  cd hotstart
  ln -sf ../fort.14
  ln -sf ../fort.13
  ln -sf ../fort.15.hotstart ./fort.15
"""
        if self._driver.wind_forcing is not None:
            if self._driver.wind_forcing.NWS in [19, 20]:
                f += "  ln -sf ../fort.15.best_track ./fort.22\n"
                f += "  aswip\n"
                f += f"  mv NWS_{self._driver.wind_forcing.NWS}_fort.22 "
                f += "fort.22\n"

        if self._executable.startswith('p'):
            if isinstance(self._server_config, SlurmConfig):
                f += "  adcprep --np $SLURM_NTASKS --partmesh\n"
                f += "  adcprep --np $SLURM_NTASKS --prepall\n"
            else:
                f += f"  adcprep --np {self._nprocs} --partmesh\n"
                f += f"  adcprep --np {self._nprocs} --prepall\n"
            f += f"  {self._mpi} {self._executable} "
        else:
            f += f"  {self._executable} "

        if not isinstance(self._server_config, SlurmConfig):
            f += f"2>&1 | tee -a ../{self._logfile}"

        f += """
  rm -rf fort.68.nc
  clean_directory
  cd ..
}
"""
        return f

    @property
    def _clean_directory(self):
        return """
clean_directory() {
  rm -rf PE*
  rm -rf partmesh.txt
  rm -rf metis_graph.txt
  rm -rf fort.13
  rm -rf fort.14
  rm -rf fort.15
  rm -rf fort.16
  rm -rf fort.80
}
"""

    @property
    def _logfile(self):
        if isinstance(self._server_config, int):
            return f'{self._executable}.log'

        if isinstance(self._server_config, SlurmConfig):
            if self._server_config._log_filename is not None:
                return self._server_config._log_filename
            else:
                return f'{self._executable}.log'

    @property
    def _executable(self):

        if self._nprocs == 1:
            if self._driver.wave_forcing is not None:
                return 'adcswan'
            else:
                return 'adcirc'
        else:
            if self._driver.wave_forcing is not None:
                if self._driver.tidal_forcing is not None:
                    return 'padcswan'
                else:
                    return 'punswan'
            else:
                return 'padcirc'

    @property
    def _mpi(self):
        if isinstance(self._server_config, SlurmConfig):
            return "srun"
        else:
            return f"mpiexec -n {self._nprocs}"

    @property
    def _server_config(self):
        return self._driver._server_config

    @property
    def _nprocs(self):
        if isinstance(self._server_config, int):
            return self._server_config
        else:
            return self._server_config.nprocs
