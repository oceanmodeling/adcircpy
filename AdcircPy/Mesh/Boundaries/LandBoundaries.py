# global imports
import numpy as np
from osgeo import osr, ogr

# local imports
from AdcircPy.Mesh.Boundaries._BaseBoundary import _BaseBoundary

# unittest imports
import unittest
import os


class LandBoundaries(_BaseBoundary):

    def __init__(self, *boundaries):
        super(LandBoundaries, self).__init__(*boundaries)

    def add_boundary(self, SpatialReference, vertices,
                     LayerName='land_boundaries', **fields):
        if isinstance(SpatialReference, int):
            EPSG = SpatialReference
            SpatialReference = osr.SpatialReference()
            SpatialReference.ImportFromEPSG(EPSG)
        else:
            assert isinstance(SpatialReference, osr.SpatialReference)
        Geometry = ogr.Geometry(ogr.wkbLineString)
        Geometry.AssignSpatialReference(SpatialReference)
        vertices = np.asarray(vertices)
        for x, y in vertices:
            Geometry.AddPoint_2D(x, y)
        super(LandBoundaries, self).add_boundary(LayerName, Geometry, **fields)


class LandBoundariesTestCase(unittest.TestCase):

    def setUp(self):
        self.LandBoundaries = LandBoundaries
        self.HSOFS_OCEAN_BOUNDARY_SHAPEFILE = os.getenv(
                                            'HSOFS_LAND_BOUNDARIES_SHAPEFILE')

    def test_empty(self):
        self.LandBoundaries()

    def test_read_Shapefile(self):
        self.LandBoundaries(self.HSOFS_LAND_BOUNDARIES_SHAPEFILE)
