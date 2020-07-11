from datetime import timedelta
import os
import pathlib
# from tempfile import TemporaryDirectory

PADCIRC_DRIVER_SCRIPT_FILENAME = pathlib.Path(os.path.dirname(__file__)).resolve() / 'padcirc_driver.sh'


class ServerConfig:
    """
    This class is used for configuring the server
    """

    def __init__(
        self,
        hostname: str = None,
        nprocs: int = None,
        wdir=None,
        binaries_prefix=None,
        port=22,
        username=None,
        password=None,
        pkey=None,
        writer_procs=None,
        source_script=None,
        additional_mpi_options=None,
        keep_wdir=False
    ):
        self._hostname = hostname
        self._nprocs = nprocs
        self._wdir = wdir
        self._binaries_prefix = binaries_prefix
        self._port = port
        self._username = username
        self._password = password
        self._pkey = pkey
        self._writer_procs = writer_procs
        self._source_script = source_script
        self._additional_mpi_options = additional_mpi_options
        self._keep_wdir = keep_wdir

    def run(
        self,
        driver,
        outdir,
        overwrite=False,
        coldstart=True,
        hotstart=True
    ):
        """
        puts driver outputs on outdir using a remote server for compute
        """
        outdir = pathlib.Path(outdir)
        if not outdir.exists():
            msg = f"{outdir} exists and overwrite is not enabled."
            raise IOError(msg)
        self.ssh.exec_command(f'mkdir -p {self.wdir}')
        self._deploy_files_to_server(driver)
        self._run_coldstart(driver)
        self._cleanup_rundir('coldstart')
        self._run_hotstart(driver)
        self._cleanup_rundir('hotstart')
        self._retrieve_files(outdir)
        if not self.keep_wdir:
            self.ssh.exec_command(f'rm -rf {self.wdir}')

    def _deploy_files_to_server(self, driver):
        outdir = tempfile.TemporaryDirectory()
        driver.dump(outdir.name)
        for item in pathlib.Path(outdir.name).glob('**/*'):
            self.sftp.put(item.absolute(), f'{self.wdir}/{item.name}')

    def _run_coldstart(self, driver):
        self._run_adcprep_command('coldstart')
        self._run_padcirc_command('coldstart', driver)

    def _run_hotstart(self, driver):
        self._run_adcprep_command('hotstart')
        self._run_padcirc_command('hotstart', driver)

    def _run_adcprep_command(self, runtype):
        cmd = f'rm -rf {self.wdir}/{runtype}; '
        cmd += f'mkdir -p {self.wdir}/{runtype}; '
        cmd += f"cd {self.wdir}/{runtype}; "
        cmd += f'ln -sf ../fort.14; '
        cmd += f'ln -sf ../fort.13; '
        cmd += f'ln -sf ../fort.15.{runtype} ./fort.15; '
        if runtype == 'hotstart':
            cmd += f'ln -sf ../coldstart/fort.67.nc ./fort.67.nc; '
        # if self.libraries_path:
        #     cmd += f'export LD_LIBRARY_PATH={self.libraries_path}:'
        #     cmd += "$LD_LIBRARY_PATH && "
        if self.source_script:
            cmd += f'source {self.source_script} && '
        cmd += f'{self.adcprep_binary} --np {self.nprocs} --partmesh && '
        cmd += f'{self.adcprep_binary} --np {self.nprocs} --prepall'
        stdin, stdout, stderr = self.ssh.exec_command(cmd)
        while True:
            out = stdout.readline()
            if not out:
                break
            print(out, end='')
        lines = stderr.readlines()
        if len(lines) > 0:
            msg = "\n"
            msg += "".join(lines)
            raise Exception(msg)

    def _run_padcirc_command(self, runtype, driver):
        cmd = ''
        if self.nprocs > 1:
            # if self.libraries_path:
            #     cmd += f'export LD_LIBRARY_PATH={self.libraries_path}:'
            #     cmd += "$LD_LIBRARY_PATH && "
            if self.source_script:
                cmd += f'source {self.source_script} && '
            cmd += f'mpiexec -n {self.nprocs} '
            if self.additional_mpi_options:
                mpi_opts = self.additional_mpi_options.strip("'\"")
                cmd += f'{mpi_opts} '
            cmd += f"--wdir {self.wdir}/{runtype} "
        cmd += f"{self.padcirc_binary}"
        self.logger.info(cmd)
        stdin, stdout, stderr = self.ssh.exec_command(cmd)
        while True:
            out = stdout.readline()
            if not out:
                break
            print(out, end='')
        lines = stderr.readlines()
        msg = "** ERROR: Elevation.gt.ErrorElev, ADCIRC stopping. **"
        if msg in "".join(lines):
            self.logger.warning(msg)
            driver._handle_blowup(lines)
        # filter IEEE_UNDERFLOW_FLAG IEEE_DENORMAL
        msg = "Note: The following floating-point exceptions are signalling:"
        lines = [line for line in lines if msg not in line]
        if len(lines) > 0:
            if msg not in "".join(lines):
                msg = "\n"
                msg += "".join(lines)
                raise Exception(msg)
            else:
                raise Exception(msg)

    def _retrieve_files(self, outdir):
        self.sftp.chdir(str(self.wdir))
        for i, walker in enumerate(self._sftp_walk(str(self.wdir))):
            if i == 0:
                parent = walker[0]
            rdir = pathlib.Path(walker[0])
            for file in walker[1]:
                rfile = rdir / file
                rsubdir = str(rdir).split(parent)[1].strip('/')
                ldir = outdir / rsubdir
                if not ldir.exists():
                    ldir.mkdir()
                self.sftp.get(
                    str(rfile),
                    str(ldir / file))

    def _sftp_walk(self, remotepath):
        """
        https://techtalkontv.wordpress.com/2016/11/05/python-pramiko-sftp-copydownload-all-files-in-a-folder-recursively-from-remote-server/
        """
        path = remotepath
        files = []
        folders = []
        for f in self.sftp.listdir_attr(remotepath):
            if S_ISDIR(f.st_mode):
                folders.append(f.filename)
            else:
                files.append(f.filename)
        if files:
            yield path, files
        for folder in folders:
            new_path = os.path.join(remotepath, folder)
            for x in self._sftp_walk(new_path):
                yield x

    def _cleanup_rundir(self, runtype):
        self.ssh.exec_command(f'rm -rf {self.wdir}/{runtype}/PE*')
        self.ssh.exec_command(f'rm -rf {self.wdir}/{runtype}/partmesh.txt')
        self.ssh.exec_command(f'rm -rf {self.wdir}/{runtype}/fort.13')
        self.ssh.exec_command(f'rm -rf {self.wdir}/{runtype}/fort.14')
        self.ssh.exec_command(f'rm -rf {self.wdir}/{runtype}/fort.15')
        if runtype == 'coldstart':
            self.ssh.exec_command(f'rm -rf {self.wdir}/{runtype}/fort.68.nc')

    @property
    def hostname(self):
        return self._hostname

    @property
    def nprocs(self):
        return self._nprocs

    @property
    def wdir(self):
        return self._wdir

    @property
    def binaries_prefix(self):
        return self._binaries_prefix

    @property
    def port(self):
        return self._port

    @property
    def username(self):
        return self._username

    @property
    def password(self):
        return self._password

    @property
    def pkey(self):
        return self._pkey

    @property
    def writer_procs(self):
        return self._writer_procs

    @property
    def source_script(self):
        return self._source_script

    @property
    def additional_mpi_options(self):
        return self._additional_mpi_options

    @property
    def keep_wdir(self):
        return self._keep_wdir

    @property
    def logger(self):
        try:
            return self.__logger
        except AttributeError:
            self.__logger = logging.getLogger(
                __name__ + '.' + self.__class__.__name__)
            return self.__logger

    @property
    def ssh(self):
        try:
            return self.__ssh
        except AttributeError:
            ssh = paramiko.SSHClient()
            ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            kwargs = {
                'hostname': self.hostname,
                'port': self.port,
                'username': self.username,
                'password': self.password
            }
            if self.pkey:
                kwargs.update({
                    'pkey': paramiko.RSAKey.from_private_key_file(
                        self.pkey)})
            # try:
            ssh.connect(**kwargs)
            # except paramiko.ssh_exception.SSHException:
            #     def auth_handler(title, _, fields):
            #         if len(fields) > 1:
            #             raise paramikoSSHException("Expecting one field only.")
            #         return [password]

            #     transport = ssh.get_transport()
            #     transport.auth_interactive(
            #         self.username,
            #         auth_handler)
            self.__ssh = ssh
            return self.__ssh

    @property
    def sftp(self):
        try:
            return self.__sftp
        except AttributeError:
            self.__sftp = self.ssh.open_sftp()
            return self.__sftp

    @property
    def padcirc_binary(self):
        if self.binaries_prefix:
            return self.binaries_prefix.absolute() / 'padcirc'
        else:
            return 'padcirc'

    @property
    def adcprep_binary(self):
        if self.binaries_prefix:
            return self.binaries_prefix.absolute() / 'adcprep'
        else:
            return 'adcprep'

    @property
    def _hostname(self):
        return self.__hostname

    @property
    def _nprocs(self):
        return self.__nprocs

    @property
    def _wdir(self):
        return self.__wdir

    @property
    def _binaries_prefix(self):
        return self.__binaries_prefix

    @property
    def _port(self):
        return self.__port

    @property
    def _username(self):
        return self.__username

    @property
    def _password(self):
        return self.__password

    @property
    def _pkey(self):
        return self.__pkey

    @property
    def _writer_procs(self):
        return self.__writer_procs

    @property
    def _source_script(self):
        return self.__source_script

    @property
    def _additional_mpi_options(self):
        return self.__additional_mpi_options

    @property
    def _keep_wdir(self):
        return self.__keep_wdir

    @_hostname.setter
    def _hostname(self, hostname):
        self.__hostname = hostname

    @_nprocs.setter
    def _nprocs(self, nprocs):
        self.__nprocs = nprocs

    @_wdir.setter
    def _wdir(self, wdir):
        if not wdir:
            wdir = f'/tmp/{uuid.uuid4().hex[:8]}'
        self.__wdir = pathlib.Path(wdir)

    @_binaries_prefix.setter
    def _binaries_prefix(self, binaries_prefix):
        if binaries_prefix:
            binaries_prefix = pathlib.Path(binaries_prefix)
        self.__binaries_prefix = binaries_prefix

    @_port.setter
    def _port(self, port):
        self.__port = port

    @_username.setter
    def _username(self, username):
        self.__username = username

    @_password.setter
    def _password(self, password):
        self.__password = password

    @_pkey.setter
    def _pkey(self, pkey):
        self.__pkey = pkey

    @_writer_procs.setter
    def _writer_procs(self, writer_procs):
        self.__writer_procs = writer_procs

    @_source_script.setter
    def _source_script(self, source_script):
        self.__source_script = source_script

    @_additional_mpi_options.setter
    def _additional_mpi_options(self, additional_mpi_options):
        self.__additional_mpi_options = additional_mpi_options

    @_keep_wdir.setter
    def _keep_wdir(self, keep_wdir):
        self.__keep_wdir = keep_wdir


# class ServerConfig:
#     def __init__(self, script: str):
#         self.script = script

#     def run(self):
#         """
#         Run the current shell script from a temporary file.
#         """

#         with TemporaryDirectory() as temporary_directory:
#             temporary_filename = os.path.join(temporary_directory, 'temp.job')
#             with open(temporary_filename, 'w') as temporary_file:
#                 temporary_file.write(self.script)
#             os.system(temporary_filename)

#     def write(self, filename: str):
#         """
#         Write shell script to the given filename.

#         :param filename: file path to shell script
#         """

#         with open(filename, 'w') as output_file:
#             output_file.write(self.script)


# class SlurmScript:  # (ServerConfig):
#     """
#     Object instance of a Slurm shell script (`*.job`).
#     """

#     def __init__(
#             self,
#             account: str,
#             slurm_ntasks: int,
#             run_name: str,
#             partition: str,
#             duration: timedelta,
#             driver_script_filename: str = None,
#             run_directory: str = '.',
#             mail_type: str = None,
#             mail_user: str = None,
#             log_filename: str = None,
#             modules: [str] = None,
#             path_prefix: str = None,
#             extra_commands: [str] = None,
#             **kwargs
#     ):
#         """
#         Instantiate a new Slurm shell script (`*.job`).

#         :param account: Slurm account name
#         :param slurm_ntasks: number of Slurm tasks
#         :param run_name: Slurm run name
#         :param partition: partition to run on
#         :param duration: time delta
#         :param driver_script_filename: file path to the driver shell script
#         :param run_directory: directory to run in
#         :param mail_type: email type
#         :param mail_user: email address
#         :param log_filename: file path to output log file
#         :param modules: list of file paths to modules to load
#         :param path_prefix: file path to prepend to the PATH
#         :param extra_commands: list of extra shell commands to insert into script
#         """

#         if driver_script_filename is None:
#             driver_script_filename = PADCIRC_DRIVER_SCRIPT_FILENAME

#         argument_sbatch_translations = {
#             'run_directory': 'D',
#             'run_name'     : 'J',
#             'account'      : 'A',
#             'log_filename' : 'output',
#             'slurm_ntasks' : 'n',
#             'duration'     : 'time',
#             'partition'    : 'partition'
#         }

#         for sbatch_argument in argument_sbatch_translations.values():
#             if sbatch_argument in kwargs:
#                 del kwargs[sbatch_argument]

#         if log_filename is None:
#             log_filename = 'sbatch.log'

#         hours, remainder = divmod(duration, timedelta(hours=1))
#         minutes, remainder = divmod(remainder, timedelta(minutes=1))
#         seconds = round(remainder / timedelta(seconds=1))
#         duration = f'{hours:02}:{minutes:02}:{seconds:02}'

#         script_prefix_lines = [
#             '#!/bin/bash --login',
#             f'#SBATCH -D {run_directory}',
#             f'#SBATCH -J {run_name}',
#             f'#SBATCH -A {account}'
#         ]

#         if mail_type is not None:
#             script_prefix_lines.append(f'#SBATCH --mail-type={mail_type}')
#         if mail_user is not None:
#             script_prefix_lines.append(f'#SBATCH --mail-user={mail_user}')
#         if log_filename is not None:
#             script_prefix_lines.append(f'#SBATCH --output={log_filename}')

#         script_prefix_lines.extend([
#             f'#SBATCH -n {slurm_ntasks}',
#             f'#SBATCH --time={duration}',
#             f'#SBATCH --partition={partition}'
#         ])

#         # append any additional SBATCH keywords and values passed to the function
#         for keyword, value in kwargs.items():
#             script_prefix_lines.append(f'#SBATCH {"-" if len(keyword) == 1 else "--"}{keyword}={value}')

#         if modules is not None:
#             for module in modules:
#                 script_prefix_lines.append(f'module load {module}')

#         if path_prefix is not None:
#             script_prefix_lines.append(f'PATH={path_prefix}:$PATH')

#         if extra_commands is not None:
#             script_prefix_lines.extend(extra_commands)

#         script_prefix_lines.append('set -e')

#         with open(driver_script_filename) as driver_script_file:
#             driver_script = ''.join([line for line in driver_script_file.readlines()][2:])

#         super().__init__('\n'.join(script_prefix_lines) + '\n\n' + driver_script)
