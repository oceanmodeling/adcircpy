#!/usr/bin/env python
import os
from datetime import datetime, timedelta
from adcircpy.mesh import AdcircMesh
from adcircpy.model import AdcircRun, TidalForcing
import unittest


class ShinnecockInletTestCase(unittest.TestCase):

    def setUp(self):
        self.mesh = AdcircMesh.open(os.getenv('SHINNECOCK_INLET_MESH'))

    def _test_planar_straight_line_graph(self):
        pslg = self.mesh.get_planar_straight_line_graph()
        pslg.make_plot()

    def _test(self):
        AdcircRun(self.mesh)

    def _test_AdcircRun_TidalForcing(self):
        end_date = datetime.now()
        start_date = end_date - timedelta(days=30.)
        _TidalForcing = TidalForcing(start_date, end_date)
        AdcircRun(self.mesh, _TidalForcing)

    def _test_set_SpatialReference(self):
        print('case 1')
        print(self.mesh.SpatialReference)
        self.mesh.make_plot(show=True)
        print('case 2')
        self.mesh._SpatialReference = 4326
        print(self.mesh.SpatialReference)
        self.mesh.make_plot(show=True)
        print('case 3')
        self.mesh.SpatialReference = 3395
        print(self.mesh.SpatialReference)
        self.mesh.make_plot(show=True)
        print('case 4')
        self.mesh._SpatialReference = None
        print(self.mesh.SpatialReference)
        self.mesh.make_plot(show=True)
        print('case 5')
        self.mesh.SpatialReference = 3395
        print(self.mesh.SpatialReference)
        self.mesh.make_plot(show=True)


if __name__ == '__main__':
    unittest.main()
