from __future__ import absolute_import, division, print_function
import os
import subprocess
import shutil
import numpy as np

def convert(xyz, ihorz, ivert,  ohorz,  overt,  vdatumdir, verbose=True):
    if os.path.isfile(vdatumdir+"/vdatum.jar")==False:
        raise FileNotFoundError("vdatum.jar not found in provided path: '{}'".format(vdatumdir))

    try: os.makedirs('vdatum_tmp/input')
    except: pass
    
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
    Popen_list.append(vdatumdir+"/vdatum.jar")
    Popen_list.append(ihorz)
    Popen_list.append(ivert)
    Popen_list.append(ohorz)
    Popen_list.append(overt)
    Popen_list.append('-nodata')
    Popen_list.append("-file:txt:space,0,1,2:vdatum_tmp/input/vdatum.txt:vdatum_tmp/output")
    # Popen_list.append("-file:txt:space,0,1,2:{}/vdatum_tmp/input/vdatum.txt:{}/vdatum_tmp/output".format(os.getcwd(),os.getcwd()))

    if verbose==True:
        print(' '.join(Popen_list))
    
    p = subprocess.Popen(Popen_list, stdout=subprocess.PIPE, bufsize=1)
    if verbose==True:
        for line in iter(p.stdout.readline, b''):
            print(line.decode("utf-8").strip('\n'),)
    p.terminate()

    xyz = np.loadtxt('vdatum_tmp/output/vdatum.txt')
    # shutil.rmtree('vdatum_tmp')
    return xyz