from AdcircPy.Mesh   import AdcircMesh
from AdcircPy.Outputs import _OutputFactory as OutputFactory

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
    def read_output(path, fort14=None, fort15=None, datum='LMSL', epsg=None, datum_grid=None):
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
        return OutputFactory(path, fort14, fort15, datum, epsg, datum_grid)

from AdcircPy.Tides import TidalDB; TidalDB.init_TPXO_cache()