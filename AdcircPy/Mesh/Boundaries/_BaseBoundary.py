# global imports
from collections.abc import Mapping
from osgeo import ogr

# unittest imports
import unittest


class _BaseBoundary(Mapping):
    __storage = dict()

    def __init__(self, *DataSource):
        for _DataSource in DataSource:
            if isinstance(_DataSource, str):
                _DataSource = ogr.Open(_DataSource)
            elif isinstance(DataSource, ogr.DataSource):
                _DataSource = DataSource
            else:
                raise Exception('Input must be an OGRString or OGRDataSource')
            for Layer in _DataSource:
                LayerName = Layer.GetName()
                if LayerName not in self._storage.keys():
                    self._storage[LayerName] = list()
                for Feature in Layer:
                    data = dict()
                    FeatureDefn = Feature.GetDefnRef()
                    if FeatureDefn.GetGeomFieldCount() != 1:
                        raise RuntimeError('Feature table must contain '
                                           + 'exactly one geometry column.')
                    data['Geometry'] = Feature.GetGeometryRef()
                    data['Fields'] = dict()
                    for i in range(Feature.GetFieldCount()):
                        FieldDefn = Feature.GetFieldDefnRef(i)
                        FieldName = FieldDefn.GetName()
                        data['Fields'][FieldName] = Feature.GetField(i)
                    self._storage[LayerName].append(data)

    def add_boundary(self, LayerName, Geometry, **fields):
        if LayerName not in self._storage.keys():
            self._storage[LayerName] = list()
        data = {'Geometry': Geometry,
                'Fields': fields}
        self._storage[LayerName].append(data)


class BaseBoundaryTestCase(unittest.TestCase):

    def setUp(self):
        self._BaseBoundary = _BaseBoundary

    def test_empty(self):
        self._BaseBoundary()
