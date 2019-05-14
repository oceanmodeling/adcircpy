import numpy as np
import getpass
# import paramiko

class ServerConfiguration(object):
    """
    This class is used for configuring the server 
    """
    def __init__(self, AdcircRun, run_directory, ADCIRC_BINARIES_PATH, numprocs=None,
                                                                hostname='localhost', port=22, username=None,                             
                                                                password=None, IdentityFile=None, hostfile=None,
                                                                module_list=None, PBS=None,
                                                                FORT14_PATH=None, FORT13_PATH=None,
                                                                WRITER_PROCS=1, environment_file=None):
        super(ServerConfiguration, self).__init__()
        self._AdcircRun = AdcircRun
        self._run_directory = run_directory
        self._ADCIRC_BINARIES_PATH = ADCIRC_BINARIES_PATH
        self._PBS = PBS
        self._hostfile = hostfile
        self._module_list = module_list
        self._hostname = hostname
        self._username = username
        self._port = port
        self._IdentityFile = IdentityFile
        self._numprocs = numprocs
        self._FORT14_PATH = FORT14_PATH
        self._FORT13_PATH = FORT13_PATH
        self._WRITER_PROCS = WRITER_PROCS
        self._environment_file = environment_file
        self.__init_numprocs()
        self.__init_FORT14_PATH()
        self.__init_FORT13_PATH()
        self.__init_ServerConfiguration()
        
    def deploy(run=True, monitor=True):
        self.__init_ssh_client()
        # self._ssh_client.

    def dump(self, path, filename='ADCIRC_RUN.sh'):
        with open(path+'/'+filename, 'w') as f:
            for line in self._ServerConfiguration:
                f.write(line)

    def printf(self):
        for line in self._ServerConfiguration:
            print(line, end="")

    def __init_numprocs(self):
        if self._PBS is not None:
            self._numprocs = self._PBS._numprocs

    def __init_FORT14_PATH(self):
        if self._FORT14_PATH is None:
            self._FORT14_PATH='{}/fort.14'.format(self._run_directory)

    def __init_FORT13_PATH(self):
        if self._FORT13_PATH is None:
            self._FORT13_PATH='{}/fort.13'.format(self._run_directory)

    def __init_ssh_client(self):
        self._ssh_client=paramiko.SSHClient()
        self._ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        try:
            self._ssh_client.connect(hostname=self._hostname,
                                                             username=self._username,
                                                             password=self._password,
                                                             key_filename=self._IdentityFile)
        except paramiko.ssh_exception.SSHException as e:
            self._transport = self._ssh_client.get_transport()
            self._transport.auth_interactive(self._username, self.__interactive_auth_handler)

    def __write_adcirc_run(self):
        self._ServerConfiguration.append('$ADCIRC_BINARIES_PATH/adcprep --np $RUN_PROCS --partmesh\n')
        self._ServerConfiguration.append('$ADCIRC_BINARIES_PATH/adcprep --np $RUN_PROCS --prepall\n')
        string='mpirun -np $NPROCS '
        if self._hostfile is not None:
            string+='--hostfile {} '.format(self._hostfile)
        if 300<=np.abs(self._AdcircRun.NWS)<=350:
            string+='$ADCIRC_BINARIES_PATH/padcswan -W $WRITER_PROCS\n'
        else:
            string+='$ADCIRC_BINARIES_PATH/padcirc -W $WRITER_PROCS\n'
        self._ServerConfiguration.append(string)
        self._ServerConfiguration.append('\n')

    def __init_ServerConfiguration(self):
        self._ServerConfiguration=list()
        if self._PBS is not None:
            for line in self._PBS._PBS:
                self._ServerConfiguration.append(line)
        else:
            self._ServerConfiguration.append("#!/bin/bash\n")
        self._ServerConfiguration.append('\n')
        self._ServerConfiguration.append('set -e\n')
        self._ServerConfiguration.append('\n')
        self._ServerConfiguration.append('# Set path to ADCIRC binaries:\n')
        self._ServerConfiguration.append('ADCIRC_BINARIES_PATH={}\n'.format(self._ADCIRC_BINARIES_PATH))
        self._ServerConfiguration.append('\n')
        self._ServerConfiguration.append('# Set PATH to the location of the mesh file and fort.13:\n')
        self._ServerConfiguration.append('FORT14_PATH={}\n'.format(self._FORT14_PATH))
        self._ServerConfiguration.append('FORT13_PATH={}\n'.format(self._FORT13_PATH))
        self._ServerConfiguration.append('\n')
        if self._module_list is not None:
            self._ServerConfiguration.append('# Load HPC modules required to run ADCIRC\n')
            for module in self._module_list:
                self._ServerConfiguration.append('module load {}\n'.format(module))
            self._ServerConfiguration.append('\n')
        if self._environment_file is not None:
            self._ServerConfiguration.append('source {}\n'.format(self._environment_file))
            self._ServerConfiguration.append('\n')
        self._ServerConfiguration.append('#--------------------------------------------#\n')
        self._ServerConfiguration.append('#-------------Begin ADCIRC run---------------#\n')
        self._ServerConfiguration.append('#--------------------------------------------#\n')
        self._ServerConfiguration.append('\n')
        self._ServerConfiguration.append('RUNDIR={}\n'.format(self._run_directory))
        self._ServerConfiguration.append('WRITER_PROCS={}\n'.format(self._WRITER_PROCS))
        if self._PBS is not None:
            self._ServerConfiguration.append('RUN_PROCS="$(($PBS_NP-$WRITER_PROCS))"\n')
        elif self._numprocs is None:
            self._ServerConfiguration.append('RUN_PROCS="$(($(nproc)-$WRITER_PROCS))"\n')
        else:
            self._ServerConfiguration.append('RUN_PROCS="$(({}-$WRITER_PROCS))"\n'.format(self._numprocs))
        self._ServerConfiguration.append('NPROCS="$(($RUN_PROCS+$WRITER_PROCS))"\n')
        self._ServerConfiguration.append('COLDSTARTDIR=$RUNDIR/coldstart\n')
        self._ServerConfiguration.append('HOTSTARTDIR=$RUNDIR/hotstart\n')
        self._ServerConfiguration.append('OUTPUTDIR=$RUNDIR/outputs\n')
        self._ServerConfiguration.append('\n')
        self._ServerConfiguration.append('#----------------run coldstart---------------#\n')
        # self._ServerConfiguration.append('rm -rf $COLDSTARTDIR\n')
        self._ServerConfiguration.append('mkdir -p $COLDSTARTDIR\n')
        self._ServerConfiguration.append('cd $COLDSTARTDIR\n')
        self._ServerConfiguration.append('ln -sf $FORT14_PATH $COLDSTARTDIR/fort.14\n')
        self._ServerConfiguration.append('ln -sf $FORT13_PATH $COLDSTARTDIR/fort.13\n')
        self._ServerConfiguration.append('ln -sf $RUNDIR/fort.15.coldstart $COLDSTARTDIR/fort.15\n')
        self.__write_adcirc_run()
        self._ServerConfiguration.append('#----------------run hotstart----------------#\n')
        # self._ServerConfiguration.append('rm -rf $HOTSTARTDIR\n')
        self._ServerConfiguration.append('mkdir -p $HOTSTARTDIR\n')
        self._ServerConfiguration.append('cd $HOTSTARTDIR\n')
        self._ServerConfiguration.append('ln -sf $COLDSTARTDIR/fort.67.nc $HOTSTARTDIR/fort.67.nc\n')
        self._ServerConfiguration.append('ln -sf $RUNDIR/fort.15.hotstart $HOTSTARTDIR/fort.15\n')
        self._ServerConfiguration.append('ln -sf $FORT14_PATH $HOTSTARTDIR/fort.14\n')
        self._ServerConfiguration.append('ln -sf $FORT13_PATH $HOTSTARTDIR/fort.13\n')
        if hasattr(self._AdcircRun, 'fort22'):
            self._ServerConfiguration.append('ln -sf $RUNDIR/fort.22.best_track $HOTSTART/fort.22.best_track\n')
            self._ServerConfiguration.append('$ADCIRC_BINARIES_PATH/aswip -n 20 -m 4 -z 2 -w fort.22.best_track\n')
            self._ServerConfiguration.append('ln -sf $HOTSTARTDIR/NWS_20_fort.22 $HOTSTARTDIR/fort.22\n')
        self.__write_adcirc_run()
        self._ServerConfiguration.append('#----------------cleanup-------------------#\n')
        # self._ServerConfiguration.append('rm -rf $OUTPUTDIR\n')
        self._ServerConfiguration.append('mkdir -p $OUTPUTDIR\n')
        self._ServerConfiguration.append('ln -f $HOTSTARTDIR/*.nc $OUTPUTDIR/\n')
        self._ServerConfiguration.append('\n')

    def __interactive_auth_handler(title, instructions, prompt_list):
        if prompt_list:
            return [getpass.getpass(prompt_list[0][0])]
        return []