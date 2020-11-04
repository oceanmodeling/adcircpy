from datetime import timedelta
import uuid

from adcircpy.server.base_config import BaseServerConfig


class SlurmConfig(BaseServerConfig):
    """
    Object instance of a Slurm shell script (`*.job`).
    """

    def __init__(
            self,
            account: str,
            ntasks: int,
            partition: str,
            walltime: timedelta,
            filename: str = 'slurm.job',
            run_directory: str = '.',
            run_name: str = None,
            mail_type: str = None,
            mail_user: str = None,
            log_filename: str = None,
            modules: [str] = None,
            path_prefix: str = None,
            extra_commands: [str] = None,
            launcher: str = 'srun',
            nodes: int = None
    ):
        """
        Instantiate a new Slurm shell script (`*.job`).

        :param account: Slurm account name
        :param ntasks: number of total tasks for Slurm to run
        :param run_name: Slurm run name
        :param partition: partition to run on
        :param walltime: time delta
        :param driver_script_filename: file path to the driver shell script
        :param run_directory: directory to run in
        :param mail_type: email type
        :param mail_user: email address
        :param log_filename: file path to output log file
        :param modules: list of file paths to modules to load
        :param path_prefix: file path to prepend to the PATH
        :param extra_commands: list of extra shell commands to insert into script
        :param launcher: command to start processes on target system (`srun`, `ibrun`, etc.)
        :param nodes: number of total nodes
        """
        self._account = account
        self._slurm_ntasks = ntasks
        self._run_name = run_name
        self._partition = partition
        self._walltime = walltime
        self._filename = filename
        self._run_directory = run_directory
        self._mail_type = mail_type
        self._mail_user = mail_user
        self._log_filename = log_filename
        self._modules = modules
        self._path_prefix = path_prefix
        self._extra_commands = extra_commands
        self._launcher = launcher
        self._nodes = nodes

    @property
    def nprocs(self):
        return self._slurm_ntasks

    @property
    def _walltime(self):
        return self.__walltime

    @_walltime.setter
    def _walltime(self, walltime):
        hours, remainder = divmod(walltime, timedelta(hours=1))
        minutes, remainder = divmod(remainder, timedelta(minutes=1))
        seconds = round(remainder / timedelta(seconds=1))
        self.__walltime = f'{hours:02}:{minutes:02}:{seconds:02}'

    @property
    def _filename(self):
        return self.__filename

    @_filename.setter
    def _filename(self, filename):
        if filename is None:
            filename = 'slurm.job'
        self.__filename = filename

    @property
    def _run_name(self):
        return self.__run_name

    @_run_name.setter
    def _run_name(self, run_name):
        if run_name is None:
            run_name = uuid.uuid4().hex
        self.__run_name = run_name

    @property
    def _run_directory(self):
        return self.__run_directory

    @_run_directory.setter
    def _run_directory(self, run_directory):
        if run_directory is None:
            run_directory = '.'
        self.__run_directory = run_directory

    @property
    def _log_filename(self):
        return self.__log_filename

    @_log_filename.setter
    def _log_filename(self, log_filename):
        if log_filename is None:
            log_filename = "slurm.log"
        self.__log_filename = log_filename

    @property
    def _prefix(self):
        f = f'#SBATCH -D {self._run_directory}\n' \
            f'#SBATCH -J {self._run_name}\n'

        if self._account is not None:
            f += f'#SBATCH -A {self._account}\n'
        if self._mail_type is not None:
            f += f'#SBATCH --mail-type={self._mail_type}\n'
        if self._mail_user is not None:
            f += f'#SBATCH --mail-user={self._mail_user}\n'
        if self._log_filename is not None:
            f += f'#SBATCH --output={self._log_filename}\n'

        f += f'#SBATCH -n {self._slurm_ntasks}\n'
        if self._nodes is not None:
            f += f'#SBATCH -N {self._nodes}\n'

        f += f'#SBATCH --time={self._walltime}\n' \
             f'#SBATCH --partition={self._partition}\n' \
             f'\n' \
             f'ulimit -s unlimited\n' \
             f'set -e\n'

        if self._modules is not None:
            f += f'\n' \
                 f'module load {" ".join(module for module in self._modules)}\n'

        if self._path_prefix is not None:
            f += f'\n' \
                 f'PATH={self._path_prefix}:$PATH\n'

        if self._extra_commands is not None:
            f += '\n'
            for command in self._extra_commands:
                f += f'{command}\n'

        return f
