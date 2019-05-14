# global imports
import numpy as np

# local imports
from AdcircPy.Mesh.Boundaries._BaseBoundary import _BaseBoundary
from AdcircPy.Mesh.UnstructuredMesh \
    import UnstructuredMesh as _UnstructuredMesh


class LandBoundaries(_BaseBoundary):

    __storage = list()
    __ibtypes = [0, 10, 20]

    def __getitem__(self, i):
        return self.storage[i]

    def __iter__(self):
        return iter(self.storage)

    def __len__(self):
        return len(self.storage)

    def get_Geometry(self, i, SpatialReference=None):
        return super(LandBoundaries, self).__get_LineStringTypeGeometry(
                                                        i, SpatialReference)

    def _add_boundary(self, UnstructuredMesh, indexes, ibtype):
        assert isinstance(UnstructuredMesh, _UnstructuredMesh)
        indexes = np.asarray(indexes)
        vertices = UnstructuredMesh.xy[indexes]
        ibtype = int(ibtype)
        assert indexes.shape[0] == vertices.shape[0]
        if ibtype not in self.ibtypes:
            raise TypeError('ibtype not valid. Allowed ibtypes are '
                            + '{}'.format(self.ibtypes))
        self.storage.append({
                        'vertices': vertices,
                        'node_indexes': indexes,
                        'SpatialReference': UnstructuredMesh.SpatialReference,
                        'ibtype': ibtype})

    @property
    def storage(self):
        return self.__storage

    @property
    def ibtypes(self):
        return self.__ibtypes
