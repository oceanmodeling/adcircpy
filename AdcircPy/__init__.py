from AdcircPy.Model import AdcircMesh
from AdcircPy.Model import ElevationStationsOutput
from AdcircPy.Tides import TPXO
from AdcircPy.Outputs import _OutputFactory

class AdcircPy(object):
  """
  Front-end class for reading ADCIRC mesh and ADCIRC outputs.
  """
  @staticmethod
  def read_mesh(fort14, fort13=None, **kwargs):
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
        datum  (string) <'MSL' (default)|'NAVD88'> : Mesh vertical datum.
    -------
    return:
    -------
        AdcirPy.AdcircMesh instance.
    """
    return AdcircMesh.from_fort14(fort14, fort13=fort13, **kwargs)
        
  @staticmethod
  def read_output(path, fort14=None, datum='LMSL', epsg=None, datum_grid=None):
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
        datum (str) <'MSL'(default)|'NAVD88'> : Mesh vertical datum.
    -------
    return:
    -------
        AdcirPy.<output>  where <output> is the output type.
    """
    return _OutputFactory(path, fort14, datum, epsg, datum_grid).get_output_instance()

