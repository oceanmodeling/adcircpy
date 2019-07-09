#!/usr/bin/env python
import os
from datetime import datetime, timedelta
from AdcircPy.mesh import AdcircMesh
from AdcircPy.model import AdcircRun, TidalForcing
import unittest


class ShinnecockInletTestCase(unittest.TestCase):

    def setUp(self):
        self._AdcircMesh = AdcircMesh.open(os.getenv('SHINNECOCK_INLET_MESH'))
        self._AdcircMesh._SpatialReference = 4326

    def test(self):
        AdcircRun(self._AdcircMesh)

    def test_AdcircRun_TidalForcing(self):
        end_date = datetime.now()
        start_date = end_date - timedelta(days=30.)
        _TidalForcing = TidalForcing(start_date, end_date)
        AdcircRun(self._AdcircMesh, _TidalForcing)


if __name__ == '__main__':
    unittest.main()
