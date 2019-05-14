# global imports
import numpy as np

# local imports
from AdcircPy.Mesh.Boundaries._BaseBoundary import _BaseBoundary
from AdcircPy.Mesh import UnstructuredMesh as _UnstructuredMesh


class OutflowBoundaries(_BaseBoundary):

    __storage = list()
    __ibtypes = [3, 13, 23]

    def __getitem__(self, i):
        return self.storage[i]

    def __iter__(self):
        return iter(self.storage)

    def __len__(self):
        return len(self.storage)

    def get_Geometry(self, i, SpatialReference=None):
        return super(OutflowBoundaries, self).__get_LineStringTypeGeometry(
                                                        i, SpatialReference)

    def _add_boundary(self, UnstructuredMesh, indexes, barrier_height,
                      supercritical_flow_coefficient, ibtype):
        assert isinstance(UnstructuredMesh, _UnstructuredMesh)
        indexes = np.asarray(indexes)
        vertices = UnstructuredMesh.xy[indexes]
        barrier_height = np.asarray(barrier_height)
        supercritical_flow_coefficient = np.asarray(
                                                supercritical_flow_coefficient)
        ibtype = int(ibtype)
        assert indexes.shape[0] == vertices.shape[0]
        assert indexes.shape[0] == barrier_height.shape[0]
        assert indexes.shape[0] == supercritical_flow_coefficient.shape[0]
        if ibtype not in self.ibtypes:
            raise TypeError('ibtype not valid. Allowed ibtypes are '
                            + '{}'.format(self.ibtypes))
        self.storage.append({
            'SpatialReference': UnstructuredMesh.SpatialReference,
            'vertices': vertices,
            'node_indexes': indexes,
            'barrier_height': barrier_height,
            'supercritical_flow_coefficient': supercritical_flow_coefficient,
            'ibtype': ibtype})

    @property
    def storage(self):
        return self.__storage

    @property
    def ibtypes(self):
        return self.__ibtypes
