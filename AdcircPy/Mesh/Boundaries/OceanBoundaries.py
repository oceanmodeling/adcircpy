# global imports
import numpy as np
from osgeo import osr, ogr

# local imports
from AdcircPy.Mesh.Boundaries._BaseBoundary import _BaseBoundary

# unittest imports
import unittest
import os


class OceanBoundaries(_BaseBoundary):

    def __init__(self, *boundaries):
        super(OceanBoundaries, self).__init__(*boundaries)

    def __getitem__(self, key):
        return self.storage['ocean_boundaries'][key]

    def __iter__(self):
        return iter(self.storage['ocean_boundaries'])

    def __len__(self):
        return len(self.storage['ocean_boundaries'])

    def add_boundary(self, SpatialReference, indexes, vertices, **fields):
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
        super(OceanBoundaries, self).add_boundary('ocean_boundaries', Geometry,
                                                  **fields)


class OceanBoundariesTestCase(unittest.TestCase):

    def setUp(self):
        self.OceanBoundaries = OceanBoundaries
        self.HSOFS_OCEAN_BOUNDARY_SHAPEFILE = os.getenv(
                                            'HSOFS_OCEAN_BOUNDARY_SHAPEFILE')

    def test_empty(self):
        self.OceanBoundaries()

    def test_read_Shapefile(self):
        self.OceanBoundaries(self.HSOFS_OCEAN_BOUNDARY_SHAPEFILE)
