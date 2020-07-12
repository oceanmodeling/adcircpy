from datetime import timedelta
import os
import pathlib
from adcircpy.model.driver import AdcircRun
from adcircpy.server import driver as adcirc_driver
# from tempfile import TemporaryDirectory


class SlurmConfig:
    """
    Object instance of a Slurm shell script (`*.job`).
    """

    def __init__(
        self,
        account: str,
        slurm_ntasks: int,
        run_name: str,
        partition: str,
        duration: timedelta,
        filename: str = 'slurm.job',
        run_directory: str = '.',
        mail_type: str = None,
        mail_user: str = None,
        log_filename: str = 'sbatch.log',
        modules: [str] = None,
        path_prefix: str = None,
        extra_commands: [str] = None,
    ):
        """
        Instantiate a new Slurm shell script (`*.job`).

        :param account: Slurm account name
        :param slurm_ntasks: number of Slurm tasks
        :param run_name: Slurm run name
        :param partition: partition to run on
        :param duration: time delta
        :param driver_script_filename: file path to the driver shell script
        :param run_directory: directory to run in
        :param mail_type: email type
        :param mail_user: email address
        :param log_filename: file path to output log file
        :param modules: list of file paths to modules to load
        :param path_prefix: file path to prepend to the PATH
        :param extra_commands: list of extra shell commands to insert into script
        """
        self._account = account
        self._slurm_ntasks = slurm_ntasks
        self._run_name = run_name
        self._partition = partition
        self._duration = duration
        self._filename = filename
        self._run_directory = run_directory
        self._mail_type = mail_type
        self._mail_user = mail_user
        self._log_filename = log_filename
        self._modules = modules
        self._path_prefix = path_prefix
        self._extra_commands = extra_commands

    def __call__(self, driver: AdcircRun):
        return f"{self._script_prefix}\n\n{adcirc_driver.bash(driver)}"

    def write(self, driver: AdcircRun, path: Union[str, pathlib]):
        with open(path, 'w') as output_file:
            output_file.write(self(driver))

    def deploy(self):
        raise NotImplementedError

    def run(self):
        raise NotImplementedError

    @property
    def _duration(self):
        return self.__duration

    @_duration.setter
    def _duration(self, duration):
        hours, remainder = divmod(duration, timedelta(hours=1))
        minutes, remainder = divmod(remainder, timedelta(minutes=1))
        seconds = round(remainder / timedelta(seconds=1))
        self.__duration = f'{hours:02}:{minutes:02}:{seconds:02}'

    @property
    def _script_prefix(self):
        script_prefix = [
            '#!/bin/bash --login',
            f'#SBATCH -D {self._run_directory}',
            f'#SBATCH -J {self._run_name}',
            f'#SBATCH -A {self._account}'
        ]
        if self._mail_type is not None:
            script_prefix.append(f'#SBATCH --mail-type={self._mail_type}')
        if self._mail_user is not None:
            script_prefix.append(f'#SBATCH --mail-user={self._mail_user}')
        if self._log_filename is not None:
            script_prefix.append(f'#SBATCH --output={self._log_filename}')

        script_prefix.extend([
            f'#SBATCH -n {self._slurm_ntasks}',
            f'#SBATCH --time={self._duration}',
            f'#SBATCH --partition={self._partition}'
        ])

        if self._modules is not None:
            for module in self._modules:
                script_prefix.append(f'module load {module}')

        if self._path_prefix is not None:
            script_prefix.append(f'PATH={self._path_prefix}:$PATH')

        if self._extra_commands is not None:
            script_prefix.extend(self._extra_commands)

        script_prefix.append('set -e')

        return '\n'.join(script_prefix)
