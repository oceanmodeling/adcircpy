# global imports
import numpy as np
from scipy.interpolate import griddata

# local imports
from AdcircPy.Datum._VDatumWrapper import _VDatumWrapper
from AdcircPy.Model import AdcircMesh as _AdcircMesh

# unittest imports
import unittest
import os
from AdcircPy.Mesh import UnstructuredMesh


class DatumMesh(UnstructuredMesh):
    """
    This class should be used to instantiate UnstructuredGrids that represents
    Datum conversion fields. Uses nearest to pad -99999.0 values. Use with
    caution.
    """

    def __init__(self, xy, values, elements, source_hdatum,
                 source_vdatum, target_hdatum, target_vdatum,
                 **kwargs):
        super(UnstructuredMesh, self).__init__(**kwargs)
        # self.xy = xy
        # self.y = y
        # self.values = values
        # self.elements = elements
        # self.source_hdatum = source_hdatum
        # self.source_vdatum = source_vdatum
        # self.target_hdatum = target_hdatum
        # self.target_vdatum = target_vdatum
        # self.description = description

    @classmethod
    def build_datum_mesh(cls, AdcircMesh, source_hdatum, target_hdatum,
                         source_vdatum, target_vdatum, vdatum_jar_path):
        if isinstance(AdcircMesh, _AdcircMesh) is False:
            AdcircMesh = _AdcircMesh.from_fort14(AdcircMesh,
                                                 vertical_datum=source_vdatum)
        xyz = _VDatumWrapper.jar_wrapper(AdcircMesh.get_xyz(), source_hdatum,
                                         source_vdatum, target_hdatum,
                                         target_vdatum, vdatum_jar_path)
        x = xyz[:, 0]
        y = xyz[:, 1]
        # VDatum uses six-digit while ADCIRC uses five digits.
        values = np.ma.masked_equal(xyz[:, 2], -999999.)
        values = AdcircMesh.z - values
        masked = np.where(values.mask)[0]
        not_masked = np.where(~values.mask)[0]
        pad_values = griddata((x[not_masked], y[not_masked]),
                              values[not_masked], (x[masked], y[masked]),
                              method='nearest')
        for i, _idx in enumerate(masked):
            values[_idx] = pad_values[i]
        xy = np.c_[x, y]
        return cls(xy, values, AdcircMesh.elements, source_hdatum,
                   source_vdatum, target_hdatum, target_vdatum,
                   AdcircMesh.description)


class DatumMeshTestCase(unittest.TestCase):

    def test_generate_DatumMesh(self):
        grid = DatumMesh.build_datum_mesh(
                    os.getenv('FORT14'),
                    source_hdatum='NAD83',
                    target_hdatum='NAD83',
                    source_vdatum="LMSL",
                    target_vdatum="MHHW",
                    vdatum_jar_path=os.getenv('VDATUM_JAR_PATH'))
        grid.dump('./fort14_lmsl2mhhw.grd')
        data = _AdcircMesh.parse_fort14('./fort14_lmsl2mhhw.grd')
        xy = np.vstack([data.pop('x'), data.pop('y')]).T
        mesh = UnstructuredMesh(xy, **data)
        mesh.make_plot(show=True)
        os.remove('./fort14_lmsl2mhhw.grd')
