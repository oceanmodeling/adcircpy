# global imports
import numpy as np
from osgeo import osr, ogr

# local imports
from AdcircPy.Mesh.Boundaries._BaseBoundary import _BaseBoundary

# unittest imports
import unittest


class InflowBoundaries(_BaseBoundary):

    def __init__(self, *boundaries):
        super(InflowBoundaries, self).__init__(*boundaries)

    def add_boundary(self, SpatialReference, vertices,
                     LayerName='inflow_boundaries', **fields):
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
        super(InflowBoundaries, self).add_boundary(LayerName, Geometry,
                                                   **fields)


class InflowBoundariesTestCase(unittest.TestCase):

    def test_empty(self):
        InflowBoundaries()
