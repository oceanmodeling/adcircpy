import os
import subprocess
import shutil
import numpy as np
import requests
import json


class VDatum(object):
  """
  This class uses lazy instantiation to retrieve the value of
  the vertical datum transformation from the VDatum REST API.
  The REST API seems to be unable to handle multiple points,
  and that's why it has been implemented using lazy instantiation.
  Potentially the return values colud be collected in the class instead.
  Very slow for large datasets.

  Example usage:
  value = VDatum(lon: float, lat: float, height: float, ihorz: str, ivert: str, ivert_units: str, overt: str, ivert_geoid='geoid12b', overt_units='m')

  Returns a float.        
  """
  
  def __new__(cls, lon: float, lat: float, height: float, ihorz: str, ivert: str, ivert_units: str, overt: str, ivert_geoid='geoid12b', overt_units='m') -> float:
    cls.__init__(cls, lon, lat, height, ihorz, ivert, ivert_units, overt, ivert_geoid, overt_units)
    return cls._height

  def __init__(self, lon, lat, height, ihorz, ivert, ivert_units, overt, ivert_geoid='geoid12b', overt_units='m'):
    """ """
    self._request(lon, lat, height, ihorz, ivert, ivert_units, overt, ivert_geoid, overt_units)

  @classmethod
  def _request(cls, lon, lat, height, ihorz, ivert, ivert_units, overt, ivert_geoid='geoid12b', overt_units='m'):
    # Reference: https://vdatum.noaa.gov/docs/services.html
    cls._url = 'https://vdatum.noaa.gov/vdatumweb/api/tidal?'
    cls._params = dict()
    cls._params['lon']       = lon
    cls._params['lat']       = lat
    cls._params['height']    = height
    cls._params['s_h_frame'] = ihorz
    cls._params['s_v_frame'] = ivert
    cls._params['s_v_unit']  = ivert_units
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
