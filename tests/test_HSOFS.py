#!/usr/bin/env python
import os
from adcircpy.mesh import AdcircMesh
import unittest


class HsofsTestCase(unittest.TestCase):

    def setUp(self):
        self.mesh = AdcircMesh.open(os.getenv('HSOFS_MESH'))
        self.mesh._SpatialReference = 4326

    def test_get_planar_straight_line_graph(self):
        pslg = self.mesh.get_planar_straight_line_graph()
        pslg.make_plot()

    def _test_ocean_boundary(self):
        self.mesh.ocean_boundaries


if __name__ == '__main__':
    unittest.main()
