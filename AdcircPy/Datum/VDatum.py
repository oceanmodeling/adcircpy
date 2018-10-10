import os
import subprocess
import shutil
import numpy as np
from AdcircPy.Tides import TidalForcing as _TidalForcing


def jar_wrapper(xyz, ihorz, ivert,  ohorz,  overt,  vdatum_jar_path, verbose=True, return_nodata=False):
    os.makedirs('vdatum_tmp/input', exist_ok=True)
    np.savetxt('vdatum_tmp/input/vdatum.txt', xyz, fmt='%.18f')
    ihorz = 'ihorz:'+ihorz
    ivert = 'ivert:'+ivert
    ohorz = 'ohorz:'+ohorz
    overt = 'overt:'+overt
    Popen_list = list()
    # Popen_list.append('export DISPLAY=:0.0')
    Popen_list.append('java')
    # Popen_list.append('-Djava.awt.headless=true')
    Popen_list.append('-jar')
    Popen_list.append(vdatum_jar_path)
    Popen_list.append(ihorz)
    Popen_list.append(ivert)
    Popen_list.append(ohorz)
    Popen_list.append(overt)
    if return_nodata==False:
        Popen_list.append('-nodata')
    Popen_list.append("-file:txt:space,0,1,2:vdatum_tmp/input/vdatum.txt:vdatum_tmp/output")
    # Popen_list.append("-file:txt:space,0,1,2:{}/vdatum_tmp/input/vdatum.txt:{}/vdatum_tmp/output".format(os.getcwd(),os.getcwd()))
    if verbose==True:
        print(' '.join(Popen_list))
    p = subprocess.Popen(Popen_list, stdout=subprocess.PIPE, bufsize=1)
    if verbose==True:
        for line in iter(p.stdout.readline, b''):
            print(line.decode("utf-8").strip('\n'),)
    p.send_signal(9)
    xyz = np.loadtxt('vdatum_tmp/output/vdatum.txt')
    shutil.rmtree('vdatum_tmp')
    return xyz
