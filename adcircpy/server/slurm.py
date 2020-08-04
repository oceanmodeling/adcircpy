from datetime import timedelta

from adcircpy.server.base_config import BaseServerConfig


class SlurmConfig(BaseServerConfig):
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
          log_filename: str = None,
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

    @property
    def nprocs(self):
        return self._slurm_ntasks

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
    def _prefix(self):
        f = f'#SBATCH -D {self._run_directory}\n' \
            f'#SBATCH -J {self._run_name}\n' \
            f'#SBATCH -A {self._account}\n'
        if self._mail_type is not None:
            f += f'#SBATCH --mail-type={self._mail_type}\n'
        if self._mail_user is not None:
            f += f'#SBATCH --mail-user={self._mail_user}\n'
        if self._log_filename is not None:
            f += f'#SBATCH --output={self._log_filename}\n'
        f += f'#SBATCH -n {self._slurm_ntasks}\n' \
             f'#SBATCH --time={self._duration}\n' \
             f'#SBATCH --partition={self._partition}\n' \
             f'\n' \
             f'set -e\n'

        if self._modules is not None:
            f += f'\n' \
                 f'module load {" ".join(module for module in self._modules)}\n'

        if self._path_prefix is not None:
            f += f'\n' \
                 f'PATH={self._path_prefix}:$PATH\n'

        if self._extra_commands is not None:
            f += f'\n'
            for command in self._extra_commands:
                f += f'{command}\n'

        return f
