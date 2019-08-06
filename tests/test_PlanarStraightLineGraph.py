#!/usr/bin/env python
import os
from adcircpy.mesh import AdcircMesh
import unittest


class PlanarStraightLineGraphTestCase(unittest.TestCase):

    def _test_hsofs(self):
        mesh = AdcircMesh.open(os.getenv('HSOFS_MESH'))
        mesh._SpatialReference = 4326
        pslg = mesh.planar_straight_line_graph
        pslg.make_plot()

    def _test_shinnecock_inlet(self):
        mesh = AdcircMesh.open(os.getenv('SHINNECOCK_INLET_MESH'))
        mesh._SpatialReference = 4326
        pslg = mesh.planar_straight_line_graph
        pslg.make_plot()

    def test_delaware_bay_oceanmesh2d(self):
        mesh = AdcircMesh.open(os.getenv('DelawareBayOceanMesh2D_MESH'))
        mesh._SpatialReference = 4326
        pslg = mesh.planar_straight_line_graph
        pslg.make_plot()


if __name__ == '__main__':
    unittest.main()
