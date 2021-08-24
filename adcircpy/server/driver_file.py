from datetime import timedelta
import os
from textwrap import indent

from adcircpy.server.base_config import BaseServerConfig
from adcircpy.server.slurm_config import SlurmConfig


class DriverFile:
    def __init__(self, driver, nprocs: int = None):
        self._driver = driver
        self.__nprocs = nprocs

    def write(self, path: str, overwrite: bool = False):
        if not os.path.exists(path) or overwrite:
            with open(path, 'w', newline='\n') as f:
                f.write(self._script)
            # os.chmod(path, 744)

    @property
    def _script(self) -> str:
        f = f'{self._shebang}\n'

        if not isinstance(self._server_config, int):
            f += self._server_config._prefix
        else:
            f += '\nulimit -s unlimited\n' '\nset -e\n'

        if self._executable.startswith('p') and isinstance(self._server_config, int):
            if self._server_config > 1:
                f += f'\nNPROCS={self._nprocs}\n'

        f += '\n'

        if self._driver.spinup_time == timedelta(seconds=0):
            f += self._single_phase_run + '\n'
        else:
            f += self._dual_phase_run + '\n'

        f += self._clean_directory + '\n' + 'main\n'

        return f

    @property
    def _shebang(self) -> str:
        f = '#!/bin/bash'
        if not isinstance(self._server_config, int):
            f += ' --login'
        return f

    @property
    def _single_phase_run(self) -> str:
        f = (
            'rm -rf work\n'
            'mkdir work\n'
            'cd work\n'
            'ln -sf ../fort.14\n'
            'ln -sf ../fort.13\n'
            'ln -sf ../fort.15 ./fort.15\n'
        )

        if self._executable.startswith('p'):
            f += (
                f'adcprep --np {self._nprocs} --partmesh\n'
                f'adcprep --np {self._nprocs} --prepall\n'
                f'{self._mpi} {self._executable}\n'
            )
        else:
            f += f'{self._executable}\n'

        f += 'clean_directory\n' + 'cd ..'

        return bash_function('main', f)

    @property
    def _dual_phase_run(self) -> str:
        return (
            self._bash_main_dual_phase
            + '\n'
            + self._run_coldstart_phase
            + '\n'
            + self._run_hotstart_phase
        )

    @property
    def _bash_main_dual_phase(self) -> str:
        error_exit_code = -1

        f = 'SECONDS=0\n' 'run_coldstart_phase\n'

        f += bash_if_statement(
            if_condition=f'grep -Rq "ERROR: Elevation.gt.ErrorElev, ADCIRC stopping." {self._logfile}',
            if_block='duration=$SECONDS\n'
            'echo "ERROR: Elevation.gt.ErrorElev, ADCIRC stopping."\n'
            'echo "Wallclock time: $(($duration / 60)) minutes and $(($duration % 60)) seconds."\n'
            f'exit {error_exit_code}',
            else_blocks=[
                'run_hotstart_phase\n'
                'duration=$SECONDS\n'
                + bash_if_statement(
                    if_condition=f'grep -Rq "ERROR: Elevation.gt.ErrorElev, ADCIRC stopping." {self._logfile}',
                    if_block='echo "ERROR: Elevation.gt.ErrorElev, ADCIRC stopping."\n'
                    'echo "Wallclock time: $(($duration / 60)) minutes and $(($duration % 60)) seconds."\n'
                    f'exit {error_exit_code}',
                ).strip('\n')
            ],
        )

        f += 'echo "Wallclock time: $(($duration / 60)) minutes and $(($duration % 60)) seconds."'

        return bash_function('main', f)

    @property
    def _run_coldstart_phase(self) -> str:
        f = (
            'rm -rf coldstart\n'
            'mkdir coldstart\n'
            'cd coldstart\n'
            'ln -sf ../fort.14\n'
            'ln -sf ../fort.13\n'
            'ln -sf ../fort.15.coldstart ./fort.15\n'
        )

        if self._executable.startswith('p'):
            if isinstance(self._server_config, SlurmConfig):
                f += (
                    'adcprep --np $SLURM_NTASKS --partmesh\n'
                    'adcprep --np $SLURM_NTASKS --prepall\n'
                )
            else:
                f += f'adcprep --np {self._nprocs} --partmesh\n'
                f += f'adcprep --np {self._nprocs} --prepall\n'
            f += f'{self._mpi} {self._executable} '
        else:
            f += f'{self._executable} '

        if not isinstance(self._server_config, SlurmConfig):
            f += f'2>&1 | tee ../{self._logfile}'

        f += '\nclean_directory\n' 'cd ..'

        return bash_function('run_coldstart_phase', f)

    @property
    def _run_hotstart_phase(self) -> str:
        f = (
            'rm -rf hotstart\n'
            'mkdir hotstart\n'
            'cd hotstart\n'
            'ln -sf ../fort.14\n'
            'ln -sf ../fort.13\n'
            'ln -sf ../fort.15.hotstart ./fort.15\n'
        )

        if self._driver.netcdf is True:
            f += 'ln -sf ../coldstart/fort.67.nc\n'
        else:
            f += 'ln -sf ../coldstart/fort.67\n'

        if self._driver.wind_forcing is not None:
            if self._driver.NWS in [17, 19, 20]:
                f += (
                    'ln -sf ../fort.22 ./fort.22\n'
                    'aswip\n'
                    f'mv NWS_{self._driver.NWS}_fort.22 fort.22\n'
                )
            else:
                msg = f'unsupported NWS value {self._driver.NWS}'
                raise NotImplementedError(msg)

        if self._executable.startswith('p'):
            if isinstance(self._server_config, SlurmConfig):
                f += (
                    'adcprep --np $SLURM_NTASKS --partmesh\n'
                    'adcprep --np $SLURM_NTASKS --prepall\n'
                )
            else:
                f += (
                    f'adcprep --np {self._nprocs} --partmesh\n'
                    f'adcprep --np {self._nprocs} --prepall\n'
                )
            f += f'{self._mpi} {self._executable} '
        else:
            f += f'{self._executable} '

        if not isinstance(self._server_config, SlurmConfig):
            f += f'2>&1 | tee -a ../{self._logfile}'
        f += '\nclean_directory\n' 'cd ..'

        return bash_function('run_hotstart_phase', f)

    @property
    def _clean_directory(self) -> str:
        return bash_function(
            'clean_directory',
            '\n'.join(
                f'rm -rf {member}'
                for member in [
                    'PE*',
                    'partmesh.txt',
                    'metis_graph.txt',
                    'fort.13',
                    'fort.14',
                    'fort.15',
                    'fort.16',
                    'fort.80',
                    'fort.68.nc',
                ]
            ),
        )

    @property
    def _logfile(self) -> str:
        if isinstance(self._server_config, int):
            return f'{self._executable}.log'

        if isinstance(self._server_config, SlurmConfig):
            if self._server_config._log_filename is not None:
                return self._server_config._log_filename
            else:
                return f'{self._executable}.log'

    @property
    def _executable(self) -> str:
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
    def _mpi(self) -> str:
        if isinstance(self._server_config, SlurmConfig):
            return self._server_config._launcher
        else:
            return f'mpiexec -n {self._nprocs}'

    @property
    def _server_config(self) -> BaseServerConfig:
        return self._driver._server_config

    @property
    def _nprocs(self) -> int:
        if self.__nprocs is not None:
            return self.__nprocs
        elif isinstance(self._server_config, int):
            return self._server_config
        else:
            return self._server_config.nprocs


def bash_if_statement(
    if_condition: str, if_block: str, else_blocks: [str] = None, indent_string: str = '  ',
) -> str:
    output = f'if {if_condition}; then\n' + indent(if_block, indent_string) + '\n'

    if else_blocks is not None:
        for else_block in else_blocks:
            if isinstance(else_block, str):
                output += 'else\n'
            else:
                assert (
                    len(else_block) == 2
                ), f'could not parse else condition / block: {else_block}'
                output += f'elif {else_block[0]}\n'
                else_block = else_block[1]
            output += indent(else_block, indent_string) + '\n'

    output += 'fi\n'

    return output


def bash_function(name: str, function_block: str, indent_string: str = '  ') -> str:
    return f'{name}() {{\n' + indent(function_block, indent_string) + '\n' + '}\n'
