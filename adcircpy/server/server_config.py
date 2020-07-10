from abc import ABC
from datetime import timedelta
import os
import pathlib
from tempfile import TemporaryDirectory

PADCIRC_DRIVER_SCRIPT_FILENAME = pathlib.Path(os.path.dirname(__file__)).resolve() / 'padcirc_driver.sh'

argument_sbatch_translations = {
    'run_directory': 'D',
    'run_name'     : 'J',
    'account'      : 'A',
    'log_filename' : 'output',
    'slurm_ntasks' : 'n',
    'duration'     : 'time',
    'partition'    : 'partition'
}


class ServerConfig(ABC):
    def __init__(self, script: str):
        self.script = script

    def run(self):
        """
        Run the current shell script from a temporary file.
        """

        with TemporaryDirectory() as temporary_directory:
            temporary_filename = os.path.join(temporary_directory, 'temp.job')
            with open(temporary_filename, 'w') as temporary_file:
                temporary_file.write(self.script)
            os.system(temporary_filename)

    def write(self, filename: str):
        """
        Write shell script to the given filename.

        :param filename: file path to shell script
        """

        with open(filename, 'w') as output_file:
            output_file.write(self.script)


class SlurmScript(ServerConfig):
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
            driver_script_filename: str = None,
            run_directory: str = '.',
            mail_type: str = None,
            mail_user: str = None,
            log_filename: str = None,
            modules: [str] = None,
            path_prefix: str = None,
            extra_commands: [str] = None,
            **kwargs
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

        if driver_script_filename is None:
            driver_script_filename = PADCIRC_DRIVER_SCRIPT_FILENAME

        for sbatch_argument in argument_sbatch_translations.values():
            if sbatch_argument in kwargs:
                del kwargs[sbatch_argument]

        if log_filename is None:
            log_filename = 'sbatch.log'

        hours, remainder = divmod(duration, timedelta(hours=1))
        minutes, remainder = divmod(remainder, timedelta(minutes=1))
        seconds = round(remainder / timedelta(seconds=1))
        duration = f'{hours:02}:{minutes:02}:{seconds:02}'

        script_prefix_lines = [
            '#!/bin/bash --login',
            f'#SBATCH -D {run_directory}',
            f'#SBATCH -J {run_name}',
            f'#SBATCH -A {account}'
        ]

        if mail_type is not None:
            script_prefix_lines.append(f'#SBATCH --mail-type={mail_type}')
        if mail_user is not None:
            script_prefix_lines.append(f'#SBATCH --mail-user={mail_user}')
        if log_filename is not None:
            script_prefix_lines.append(f'#SBATCH --output={log_filename}')

        script_prefix_lines.extend([
            f'#SBATCH -n {slurm_ntasks}',
            f'#SBATCH --time={duration}',
            f'#SBATCH --partition={partition}'
        ])

        # append any additional SBATCH keywords and values passed to the function
        for keyword, value in kwargs.items():
            script_prefix_lines.append(f'#SBATCH {"-" if len(keyword) == 1 else "--"}{keyword}={value}')

        if modules is not None:
            for module in modules:
                script_prefix_lines.append(f'module load {module}')

        if path_prefix is not None:
            script_prefix_lines.append(f'PATH={path_prefix}:$PATH')

        if extra_commands is not None:
            script_prefix_lines.extend(extra_commands)

        script_prefix_lines.append('set -e')

        with open(driver_script_filename) as driver_script_file:
            driver_script = ''.join([line for line in driver_script_file.readlines()][2:])

        super().__init__('\n'.join(script_prefix_lines) + '\n\n' + driver_script)
