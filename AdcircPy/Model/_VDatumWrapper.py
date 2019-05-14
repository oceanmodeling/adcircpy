import os
import subprocess
import shutil
import numpy as np
import requests
import json


from AdcircPy.utils import get_cache_dir


class _VDatumWrapper(object):
    """
    This class retrieves the value of the vertical datum transformation from
    the VDatum REST API.
    """

    __params = {
    
    
    }

    def __new__(cls, lon: float, lat: float, height: float, ihorz: str,
                ivert: str, ivert_units: str, overt: str,
                ivert_geoid='geoid12b', overt_units='m') -> float:

        cls.__init__(cls, lon, lat, height, ihorz, ivert, ivert_units, overt,
                     ivert_geoid, overt_units)
        return cls._height

    def __init__(self, lon, lat, height, ihorz, ivert, ivert_units, overt,
                 ivert_geoid='geoid12b', overt_units='m'):
        """ """
        self._request(lon, lat, height, ihorz, ivert, ivert_units, overt,
                      ivert_geoid, overt_units)

    @classmethod
    def request(cls, lon, lat, height, ihorz, ivert, ivert_units, overt,
                ivert_geoid='geoid12b', overt_units='m'):
        # Reference: https://vdatum.noaa.gov/docs/services.html
        cls._url = 'https://vdatum.noaa.gov/vdatumweb/api/tidal?'
        cls._params = dict()
        cls._params['lon'] = lon
        cls._params['lat'] = lat
        cls._params['height'] = height
        cls._params['s_h_frame'] = ihorz
        cls._params['s_v_frame'] = ivert
        cls._params['s_v_unit'] = ivert_units
        cls._params['s_v_geoid'] = ivert_geoid
        cls._params['t_v_frame'] = overt
        cls._params['t_v_units'] = overt_units
        response = requests.get(cls._url, params=cls._params)
        response.raise_for_status()
        json_data = json.loads(response.text)
        if 'errorCode' in json_data.keys():
            raise Exception(response.text)
        cls._height = float(json_data['tar_height'])

    @staticmethod
    def jar_wrapper(xyz, ihorz, ivert,  ohorz,  overt,  vdatum_jar_path,
                    verbose=True):
        """
        read from stdout raises "does not support long argument list"
        Writing output to hard-drive instead.
        jar input field for reading from stdout: 'pt:{},{},{}'.format(x, y, z))
        """
        assert os.path.isfile(vdatum_jar_path)
        popen = ['java']
        popen.append('-jar')
        popen.append('{}'.format(vdatum_jar_path))
        popen.append('ihorz:{}'.format(ihorz))
        popen.append('ivert:{}'.format(ivert))
        popen.append('ohorz:{}'.format(ohorz))
        popen.append('overt:{}'.format(overt))
        vdatum_tmp = get_cache_dir() + '/vdatum_tmp'
        # VDatum only outputs lower-cased paths.
        vdatum_tmp = vdatum_tmp.lower()
        os.makedirs(vdatum_tmp + '/input', exist_ok=True)
        np.savetxt(vdatum_tmp + '/input/vdatum.txt', xyz, fmt='%.18f')
        popen.append("-file:txt:space,0,1,2:"
                     + "{}/input/vdatum.txt".format(vdatum_tmp)
                     + ":{}/output".format(vdatum_tmp))
        p = subprocess.Popen(popen, stdout=subprocess.PIPE, bufsize=1)
        if verbose is True:
            print()
            for line in iter(p.stdout.readline, b''):
                print(line.decode("utf-8").strip('\n'))
        p.send_signal(9)
        xyz = np.loadtxt('{}/output/vdatum.txt'.format(vdatum_tmp))
        shutil.rmtree('{}'.format(vdatum_tmp))
        return xyz
