# global imports
import numpy as np

# local imports
from AdcircPy.Mesh.Boundaries._BaseBoundary import _BaseBoundary
from AdcircPy.Mesh.UnstructuredMesh \
    import UnstructuredMesh as _UnstructuredMesh


class CulvertBoundaries(_BaseBoundary):

    __storage = list()
    __ibtypes = [5, 25]

    def __getitem__(self, i):
        return self.storage[i]

    def __iter__(self):
        return iter(self.storage)

    def __len__(self):
        return len(self.storage)

    def get_Geometry(self, i, SpatialReference=None):
        return super(CulvertBoundaries,
                     self).__get_MultiLineStringTypeGeometry(
                                                        i, SpatialReference)

    def _add_boundary(self, UnstructuredMesh, front_face_indexes,
                      back_face_indexes, barrier_heights,
                      subcritical_flow_coefficients,
                      supercritical_flow_coefficients,
                      cross_barrier_pipe_heights, friction_factors,
                      pipe_diameters, ibtype):
        assert isinstance(UnstructuredMesh, _UnstructuredMesh)
        front_face_indexes = np.asarray(front_face_indexes)
        back_face_indexes = np.asarray(back_face_indexes)
        front_face_vertices = UnstructuredMesh.xy[front_face_indexes]
        back_face_vertices = UnstructuredMesh.xy[back_face_indexes]
        barrier_heights = np.asarray(barrier_heights)
        subcritical_flow_coefficients = np.asarray(
                                                subcritical_flow_coefficients)
        supercritical_flow_coefficients = np.asarray(
                                            supercritical_flow_coefficients)
        cross_barrier_pipe_heights = np.asarray(cross_barrier_pipe_heights)
        friction_factors = np.asarray(friction_factors)
        pipe_diameters = np.asarray(pipe_diameters)
        ibtype = int(ibtype)
        indexes = front_face_indexes
        assert indexes.shape == back_face_indexes.shape
        assert front_face_indexes.shape[0] == barrier_heights.shape[0]
        assert front_face_indexes.shape[0] \
            == subcritical_flow_coefficients.shape[0]
        assert front_face_indexes.shape[0] \
            == supercritical_flow_coefficients.shape[0]
        assert front_face_indexes.shape[0]  \
            == cross_barrier_pipe_heights.shape[0]
        assert front_face_indexes.shape[0] == friction_factors.shape[0]
        assert front_face_indexes.shape[0] == pipe_diameters.shape[0]
        if ibtype not in self.ibtypes:
            raise TypeError('ibtype not valid. Allowed ibtypes are '
                            + '{}'.format(self.ibtypes))
        self.storage.append({
            'SpatialReference': UnstructuredMesh.SpatialReference,
            'front_face_vertices': front_face_vertices,
            'back_face_vertices': back_face_vertices,
            'front_face_indexes': front_face_indexes,
            'back_face_indexes': back_face_indexes,
            'barrier_heights': barrier_heights,
            'subcritical_flow_coefficients': subcritical_flow_coefficients,
            'supercritical_flow_coefficients': supercritical_flow_coefficients,
            'cross_barrier_pipe_heights': cross_barrier_pipe_heights,
            'friction_factors': friction_factors,
            'pipe_diameters': pipe_diameters,
            'ibtype': ibtype})

    @property
    def storage(self):
        return self.__storage

    @property
    def ibtypes(self):
        return self.__ibtypes
