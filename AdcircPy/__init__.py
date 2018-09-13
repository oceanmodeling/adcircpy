from AdcircPy.Mesh   import AdcircMesh
from AdcircPy.Outputs import _OutputFactory

class AdcircPy(object):

    @staticmethod
    def read_mesh(fort14, **kwargs):
        """
        Reads ADCIRC input files.
        
        -----------
        Parameters:
        -----------
            fort14 (string) :  Path to fort.14 file. 
        
        -------------------
        Optional arguments:
        -------------------
            fort13 (string) <None> : Path to optional fort.13 file.
            fort15 (string) <None> : Path to optional fort.15 file.
            datum  (string) <'MSL' (default)|'NAVD88'> : Mesh vertical datum.
        
        -------
        return:
        -------
            AdcirPy.Mesh instance.
        """
        return AdcircMesh.from_fort14(fort14, **kwargs)
        
    @staticmethod
    def read_output(path, **kwargs):
        """
        Reads ADCIRC output files. Supports both ASCII and NetCDF.
        
        -----------
        Parameters:
        -----------
            path (str): Path to output file.
        
        -------------------
        Optional arguments:
        -------------------
            fort14 (str) <None> : Path to fort.14 file.
                                  Required for ASCII gridded field outputs (e.g. maxele).
                                  Optional in case of NetCDF files (used for inclusion of boundary data).
            fort13 (str) <None> : Path to optional fort.13 file.
            fort15 (str) <None> : Path to optional fort.15 file.
            datum (str) <'MSL'(default)|'NAVD88'> : Mesh vertical datum.
        
        -------
        return:
        -------
            AdcirPy.<output>  where <output> is the output type.
        """
        return _OutputFactory(path, **kwargs)
    

import os
import wget
import tarfile
_cachedir = os.getenv('LOCALAPPDATA')
if _cachedir is None:
    _cachedir = os.getenv('HOME')+'/.cache/AdcircPy'
else: 
    _cachedir += '/AdcircPy'
os.makedirs(_cachedir, exist_ok=True)
if os.path.isfile(_cachedir+"/h_tpxo9.v1.nc")==False:
    print('Building TPXO database cache on {}, please wait...'.format(_cachedir+"/h_tpxo9.v1.nc"))
    print('(This will only happen the first time you run this software)')
    url='ftp://ftp.oce.orst.edu/dist/tides/Global/tpxo9_netcdf.tar.gz'
    if os.path.isfile(_cachedir+"/tpxo9_netcdf.tar.gz")==False:
        tpxo=wget.download(url, out=_cachedir+"/tpxo9_netcdf.tar.gz")
        tpxo=tarfile.open(tpxo)
    else:
        tpxo=tarfile.open(tpxo)
