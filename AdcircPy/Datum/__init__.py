from AdcircPy.Datum import _DatumGrid
from AdcircPy.Datum import _VDatum
import requests
import json


class DatumGrid(object):
    """
    This class should be used to instantiate UnstructuredGrids that 
    represent Datum 
    """

    def __init__(self, **kwargs):
        self.x = kwargs.pop("x", None)
        self.y = kwargs.pop("y", None)
        self.values = kwargs.pop("values", None)
        self._type = kwargs.pop("_type", None)

    @staticmethod
    def msl_to_navd88(Datum_grid):
        return _DatumGrid.msl_to_navd88(Datum_grid)

    def convert(self, Mesh, **kwargs):
        return _DatumGrid.convert(self, Mesh, **kwargs)

class VDatum(object):
    """
    This class uses lazy instantiation to retrieve the value of
    the vertical datum transformation from the VDatum REST API.
    The REST API seems to be unable to handle multiple points,
    and that's why it has been implemented using lazy instantiation.
    Potentially the return values colud be collected in the class instead.
    Very slow for large datasets.
    It is advisable to filter 

    Example usage:
        value = VDatum(lon: float, lat: float, height: float, ihorz: str, ivert: str, ivert_units: str, overt: str, ivert_geoid='geoid12b', overt_units='m')

    Returns a float.        
    """
    
    def __new__(cls, lon: float, lat: float, height: float, ihorz: str, ivert: str, ivert_units: str, overt: str, ivert_geoid='geoid12b', overt_units='m') -> float:
        cls.__init__(cls, lon, lat, height, ihorz, ivert, ivert_units, overt, ivert_geoid, overt_units)
        return cls._height

    def __init__(self, lon, lat, height, ihorz, ivert, ivert_units, overt, ivert_geoid='geoid12b', overt_units='m'):
        self._request(lon, lat, height, ihorz, ivert, ivert_units, overt, ivert_geoid, overt_units)

    @staticmethod
    def jar_wrapper(xyz, ihorz, ivert,  ohorz,  overt,  vdatum_jar_dir, **kwargs):
        return _VDatum.jar_wrapper(xyz, ihorz, ivert,  ohorz,  overt,  vdatum_jar_dir, **kwargs)

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