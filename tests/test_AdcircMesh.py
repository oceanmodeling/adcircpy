#!/usr/bin/env python
import os
from adcircpy.mesh import AdcircMesh
import unittest


class ShinnecockInletTestCase(unittest.TestCase):

    def setUp(self):
        self.mesh = AdcircMesh.open(os.getenv('SHINNECOCK_INLET_MESH'))
        self.mesh._SpatialReference = 4326

    def _test_transform_SpatialReference(self):
        self.mesh.SpatialReference = 3395
        self.mesh.make_plot(True)

    def _test_remove_ocean_boundary(self):
        print(self.mesh.ocean_boundaries)
        self.mesh.remove_ocean_boundary(0)
        print(self.mesh.ocean_boundaries)

    def _test_write_fort14(self):
        self.mesh.write_fort14()


class HsofsTestCase(unittest.TestCase):

    def setUp(self):
        self.mesh = AdcircMesh.open(os.getenv('HSOFS_MESH'))
        self.mesh._SpatialReference = 4326

    def test_transform_SpatialReference(self):
        # self.mesh.SpatialReference = 3395
        self.mesh.make_plot(show=True)

    def _test_print_boundaries(self):
        print(self.mesh.ocean_boundaries)
        print(self.mesh.land_boundaries)
        print(self.mesh.inner_boundaries)
        print(self.mesh.weir_boundaries)

    def _test_import_attributes(self):
        self.mesh.import_attributes(os.getenv('HSOFS_ATTRIBUTES'))

    def _test_write_fort14(self):
        self.mesh.write_fort14()


if __name__ == '__main__':
    unittest.main()
