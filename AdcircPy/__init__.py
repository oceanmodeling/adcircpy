from AdcircPy.Model import AdcircMesh
from AdcircPy.Model import ElevationStationsOutput
from AdcircPy.Model import ElevationGlobalOutput
from AdcircPy.Model import VelocityStationsOutput
from AdcircPy.Model import VelocityGlobalOutput
from AdcircPy.DEM import DEM
from AdcircPy.Tides import TPXO
from AdcircPy.Outputs import _OutputFactory

class AdcircPy(object):
  """
  Front-end class for reading ADCIRC mesh and ADCIRC outputs.
  """
  @staticmethod
  def read_mesh(fort14, datum=None, epsg=4326, fort13=None, datum_grid=None):
    """
    Reads ADCIRC input files.
    -----------
    Parameters:
    -----------
        fort14 (string)     : Path to fort.14 file. 
    -------------------
    Optional arguments:
    -------------------
        fort13 (string)     : Path to optional fort.13 file.
        datum  (string)     : Mesh vertical datum. Usually 'LMSL' or 'NAVD88'.
        epsg   (int)        : Corresponds to epsg projection. Use 4326 for WGS84 (default)
        datum_grid (string) : Path to datum conversion grid. Useful when doing HWM validations.
    -------
    return:
    -------
        AdcirPy.AdcircMesh instance.
    """
    return AdcircMesh.from_fort14(fort14, datum, epsg, fort13, datum_grid)
        
  @staticmethod
  def read_output(path, fort14=None, datum=None, epsg=4326, datum_grid=None, fort15=None):
    """
    Reads ADCIRC output files and returns the appropriate output type class.
    Supports ASCII and NetCDF outputs.
    -----------
    Parameters:
    -----------
        path (str)          : Path to output file.
    -------------------
    Optional arguments:
    -------------------
        fort14 (str)        : Path to fort.14 file. Required for ASCII gridded field outputs (e.g. maxele).
                              Optional in case of NetCDF files (used for inclusion of boundary data).
        datum  (string)     : Mesh vertical datum. Usually 'LMSL' or 'NAVD88'.
        epsg   (int)        : Corresponds to epsg projection. Use 4326 for WGS84 (default)
        datum_grid (string) : Path to datum conversion grid. Useful when doing HWM validations.
    -------
    return:
    -------
        AdcirPy.<output>  where <output> is the output type.
    """
    return _OutputFactory(path, fort14, datum, epsg, datum_grid, fort15).get_output_instance()

