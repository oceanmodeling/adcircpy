import numpy as np
from adcircpy.mesh.base import EuclideanMesh2D


class NodalAttributes:

    def __init__(self, mesh, fort13=None):
        self._mesh = mesh
        self._fort13 = fort13

    def _import_fort13(self, fort13):
        fort13 = self.parse_fort13(fort13)
        if fort13.pop('NumOfNodes') != len(self.node_id):
            raise Exception('fort.13 file does not match the mesh.')
        self._AGRID = fort13.pop('AGRID')
        for attribute, data in fort13.items():
            values = np.asarray(data['values'])
            if values.ndim == 1:
                values = values.reshape((values.shape[0], 1))
            full_values = np.full(
                (self.values.size,
                    np.asarray(data['default_values']).flatten().size),
                np.nan)
            for i, idx in enumerate(data['indexes']):
                for j, value in enumerate(values[i, :].tolist()):
                    full_values[idx, j] = value
            idxs = np.where(np.isnan(full_values).all(axis=1))[0]
            for idx in idxs:
                for i, value in enumerate(data['default_values']):
                    full_values[idx, i] = value
            # converts from column major to row major, leave it column major.
            # if full_values.shape[1] == 1:
            #     full_values = full_values.flatten()
            self.add_nodal_attribute(attribute, data['units'])
            self.set_nodal_attribute(attribute, full_values)

    @property
    def _mesh(self):
        return self.__mesh

    @_mesh.setter
    def _mesh(self, mesh):
        assert isinstance(mesh, EuclideanMesh2D)



    @mannings_n_at_sea_floor.setter
    def mannings_n_at_sea_floor(self, mannings_n_at_sea_floor):
        self.add_nodal_attribute('mannings_n_at_sea_floor', 'meters')
        self.set_nodal_attribute(
            'mannings_n_at_sea_floor', mannings_n_at_sea_floor)
