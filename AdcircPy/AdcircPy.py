from AdcircPy.Model import AdcircMesh
from AdcircPy.Outputs._OutputFactory import _OutputFactory


def read_mesh(fort14, SpatialReference=4326, vertical_datum=None,
              datum_mesh=None):
    """
    Reads ADCIRC input files.
    -----------
    Parameters:
    -----------
        fort14: Path to fort.14 file.

    -------------------
    Optional arguments:
    -------------------
        fort13: Path to optional fort.13 file.

        datum: Mesh vertical datum. Usually 'LMSL'

        SpatialReference (int): Corresponds to epsg projection.
    -------
    returns:
    -------
        Instance of AdcirPy.Model.AdcircMesh
    """

    return AdcircMesh(fort14, SpatialReference, vertical_datum, datum_mesh)


def read_output(path, fort14=None, SpatialReference=4326, vertical_datum=None,
                datum_mesh=None):
    """
    Reads ADCIRC output files and returns the appropriate output type class.
    Supports ASCII and NetCDF outputs.
    -----------
    Parameters:
    -----------
        path: Path to output file.
    -------------------
    Optional arguments:
    -------------------
        Same are read_mesh()
    -------
    return:
    -------
        ChildInstance of appropriate output type.
    """
    return _OutputFactory(path, fort14, SpatialReference, vertical_datum,
                          datum_mesh).output
