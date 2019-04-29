# global imports
import numpy as np
from osgeo import osr, ogr

# local imports
from AdcircPy.Mesh.Boundaries._BaseBoundary import _BaseBoundary

# unittest imports
import unittest
import os


class WeirBoundaries(_BaseBoundary):

    def __init__(self, *boundaries):
        super(WeirBoundaries, self).__init__(*boundaries)

    def add_boundary(self, SpatialReference, front_face, back_face,
                     LayerName='weir_boundaries', **fields):
        if isinstance(SpatialReference, int):
            EPSG = SpatialReference
            SpatialReference = osr.SpatialReference()
            SpatialReference.ImportFromEPSG(EPSG)
        else:
            assert isinstance(SpatialReference, osr.SpatialReference)
        FrontFaceGeometry = ogr.Geometry(ogr.wkbLineString)
        BackFaceGeometry = ogr.Geometry(ogr.wkbLineString)
        Geometry = ogr.Geometry(ogr.wkbMultiLineString)
        Geometry.AssignSpatialReference(SpatialReference)
        front_face = np.asarray(front_face)
        back_face = np.asarray(back_face)
        assert front_face.shape == back_face.shape
        for x, y in front_face:
            FrontFaceGeometry.AddPoint_2D(x, y)
        for x, y in back_face:
            BackFaceGeometry.AddPoint_2D(x, y)
        Geometry.AddGeometry(FrontFaceGeometry)
        Geometry.AddGeometry(BackFaceGeometry)
        super(WeirBoundaries, self).add_boundary(LayerName, Geometry, **fields)


class WeirBoundariesTestCase(unittest.TestCase):

    def setUp(self):
        self.HSOFS_WEIR_BOUNDARIES_SHAPEFILE = os.getenv(
                                            'HSOFS_WEIR_BOUNDARIES_SHAPEFILE')

    def test_empty(self):
        WeirBoundaries()

    def test_read_Shapefile(self):
        WeirBoundaries(self.HSOFS_WEIR_BOUNDARIES_SHAPEFILE)
