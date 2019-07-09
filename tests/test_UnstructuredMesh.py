#!/usr/bin/env python
import os
import numpy as np
from adcircpy.mesh import UnstructuredMesh, AdcircMesh
import unittest


class UnstructuredMeshTestCase(unittest.TestCase):

    def setUp(self):
        fort14 = AdcircMesh.parse_fort14(os.getenv('SHINNECOCK_INLET_MESH'))
        vertices = np.vstack([fort14['x'], fort14['y']]).T
        self.mesh = UnstructuredMesh(
            vertices, fort14['elements'], fort14['values'], 4326)

    def _test_change_SpatialReference(self):
        self.mesh.SpatialReference = 3395
        self.mesh.make_plot(True)

    def test_change_SpatialReference_PSLG(self):
        # self.mesh.SpatialReference = 3395
        self.mesh.planar_straight_line_graph.make_plot()


if __name__ == '__main__':
    unittest.main()
