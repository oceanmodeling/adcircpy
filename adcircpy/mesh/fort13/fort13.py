from adcircpy.mesh.base import EuclideanMesh2D


class NodalAttributes:

    def __init__(self, mesh):
        self._mesh = mesh

    @property
    def _mesh(self):
        return self.__mesh

    @_mesh.setter
    def _mesh(self, mesh):
        assert isinstance(mesh, EuclideanMesh2D)
