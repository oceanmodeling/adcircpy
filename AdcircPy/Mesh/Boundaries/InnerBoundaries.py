# global imports
import numpy as np
from osgeo import osr, ogr

# local imports
from AdcircPy.Mesh.Boundaries._BaseBoundary import _BaseBoundary

# unittest imports
import unittest
import os


class InnerBoundaries(_BaseBoundary):

    def __init__(self, *boundaries):
        super(InnerBoundaries, self).__init__(*boundaries)

    def add_boundary(self, SpatialReference, vertices,
                     LayerName='inner_boundaries', **fields):
        if isinstance(SpatialReference, int):
            EPSG = SpatialReference
            SpatialReference = osr.SpatialReference()
            SpatialReference.ImportFromEPSG(EPSG)
        else:
            assert isinstance(SpatialReference, osr.SpatialReference)
        Geometry = ogr.Geometry(ogr.wkbLinearRing)
        Geometry.AssignSpatialReference(SpatialReference)
        vertices = np.asarray(vertices)
        for x, y in vertices:
            Geometry.AddPoint_2D(x, y)
        super(InnerBoundaries, self).add_boundary(LayerName, Geometry,
                                                  **fields)


class InnerBoundariesTestCase(unittest.TestCase):

    def setUp(self):
        self.InnerBoundaries = InnerBoundaries
        self.HSOFS_OCEAN_BOUNDARY_SHAPEFILE = os.getenv(
                                            'HSOFS_INNER_BOUNDARIES_SHAPEFILE')

    def test_empty(self):
        self.InnerBoundaries()

    def test_read_Shapefile(self):
        self.InnerBoundaries(self.HSOFS_INNER_BOUNDARIES_SHAPEFILE)
