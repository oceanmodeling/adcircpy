#!/usr/bin/env python
import os
import numpy as np
from adcircpy.mesh import AdcircMesh
import unittest


class HsofsTestCase(unittest.TestCase):

    def setUp(self):
        self.mesh = AdcircMesh.open(os.getenv('HSOFS_MESH'))
        self._mesh._SpatialReference = 4326

    
    def test_Isabel2003(self):
        pass


    def _test_pslg(self):
        # self.mesh.SpatialReference = 3395
        import time
        start = time.time()
        pslg = self.mesh.planar_straight_line_graph
        print('{}'.format(time.time()-start))
        pslg.make_plot()


if __name__ == '__main__':
    unittest.main()
