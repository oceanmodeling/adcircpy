# global imports
import numpy as np
from osgeo import osr, ogr

# local imports
from AdcircPy.Mesh.Boundaries._BaseBoundary import _BaseBoundary

# unittest imports
import unittest


class OutflowBoundaries(_BaseBoundary):

    def __init__(self, *boundaries):
        super(OutflowBoundaries, self).__init__(*boundaries)

    def add_boundary(self, SpatialReference, vertices,
                     LayerName='outflow_boundaries', **fields):
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
        super(OutflowBoundaries, self).add_boundary(LayerName, Geometry,
                                                    **fields)


class OutflowBoundariesTestCase(unittest.TestCase):

    def test_empty(self):
        OutflowBoundaries()
