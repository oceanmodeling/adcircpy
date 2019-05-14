# global imports
import numpy as np

# local imports
from AdcircPy.Mesh.Boundaries._BaseBoundary import _BaseBoundary
from AdcircPy.Mesh.UnstructuredMesh \
    import UnstructuredMesh as _UnstructuredMesh


class OceanBoundaries(_BaseBoundary):

    __storage = list()

    def __getitem__(self, i):
        return self.storage[i]

    def __iter__(self):
        return iter(self.storage)

    def __len__(self):
        return len(self.storage)

    def get_Geometry(self, i, SpatialReference=None):
        return super(OceanBoundaries, self).__get_LineStringTypeGeometry(
                                                        i, SpatialReference)

    def _add_boundary(self, UnstructuredMesh, indexes, ibtype=None):
        assert isinstance(UnstructuredMesh, _UnstructuredMesh)
        indexes = np.asarray(indexes)
        vertices = UnstructuredMesh.xy[indexes, :]
        assert indexes.shape[0] == vertices.shape[0]
        self.storage.append({
                        'vertices': vertices,
                        'node_indexes': indexes,
                        'SpatialReference': UnstructuredMesh.SpatialReference})

    @property
    def storage(self):
        return self.__storage
