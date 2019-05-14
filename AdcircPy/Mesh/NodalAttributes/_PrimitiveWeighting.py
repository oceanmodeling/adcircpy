# global imports
from collections import defaultdict
from itertools import permutations
import numpy as np
from matplotlib.tri import Triangulation


# local imports
from AdcircPy.Mesh.NodalAttributes._BaseNodalAttribute \
    import _BaseNodalAttribute


class _PrimitiveWeighting(_BaseNodalAttribute):

    __name = 'primitive_weighting_in_continuity_equation'

    def __init__(self, UnstructuredMesh, values=None, TAU0=-3,
                 default_values=None, units='unitless', **kwargs):
        if TAU0 == 3:
            values, default_values = self.tau0_gen(UnstructuredMesh, **kwargs)
        super(_PrimitiveWeighting, self).__init__(UnstructuredMesh, values,
                                                  default_values, units)
        self._TAU0 = TAU0

    def auto_generate_indexes(self, UnstructuredMesh, default_value=0.03,
                              **kwargs):
        self.__values = self.tau0_gen(UnstructuredMesh, **kwargs)
        self.__default_values = np.asarray(default_value)
        self.__TAU0 = 3.

    def __set_Tau0FullDomainMin(self, Tau0FullDomainMin):
        self.__Tau0FullDomainMin = float(Tau0FullDomainMin)

    def __set_Tau0FullDomainMax(self, Tau0FullDomainMax):
        self.__Tau0FullDomainMax = float(Tau0FullDomainMax)

    @classmethod
    def from_TAU0(cls, UnstructuredMesh, TAU0, **kwargs):
        TAU0 = float(TAU0)
        if TAU0 >= 0. and TAU0 <= 1.:
            PrimitiveWeighting = cls(UnstructuredMesh, None, TAU0)
        elif TAU0 == 3.:
            PrimitiveWeighting = cls(UnstructuredMesh, None, TAU0)
        elif TAU0 in [-1., -2.]:
            PrimitiveWeighting = cls(UnstructuredMesh, None, TAU0)
        elif TAU0 == -3:
            PrimitiveWeighting = cls._from_fort13(UnstructuredMesh, **kwargs)
        elif TAU0 == -5.:
            PrimitiveWeighting = cls(UnstructuredMesh, None, TAU0)
            Tau0FullDomainMin = kwargs.pop('Tau0FullDomainMin', 0.005)
            PrimitiveWeighting.__set_Tau0FullDomainMin(Tau0FullDomainMin)
            Tau0FullDomainMax = kwargs.pop('Tau0FullDomainMax', 0.02)
            PrimitiveWeighting.__set_Tau0FullDomainMax(Tau0FullDomainMax)

        else:
            raise TypeError(
                'Unrecognized value for TAU0. Allowed values values are: \n'
                + '\t3: (default) Auto calculates the TAU0 values internally.'
                + '\n\t0: pure wave equation.\n'
                + '\t(0, 1): pure primitive continuity equation\n'
                + '\t-1: TAU0 calulated from depth.\n'
                + '\t-2: TAU0 calulated from depth ranges (see ADCIRC manual).'
                + '\n\t-3: Keeps the TAU0 founds in the fort.13\n'
                + '\t-5: From the ADCIRC manual: "TAU0 varies spatially and in'
                + 'time, and is dependent on the local friction" Also requires'
                + " specification of 'Tau0FullDomainMin and Tau0FullDomainMax "
                + 'in kwargs.\n For more information see: https://adcirc.org'
                + '/home/documentation/users-manual-v53/parameter-definitions/'
                + '#TAU0')
        PrimitiveWeighting._TAU0 = TAU0
        return PrimitiveWeighting

    @staticmethod
    def tau0_gen(UnstructuredMesh, default_value=0.03,
                 threshold_distance=1750., shallow_tau0=0.02, deep_tau0=0.005,
                 threshold_depth=-10.):
        """
        Reimplementation of tau0_gen.f by Robert Weaver (2008)
        1) computes  distance to each neighboring node
        2) averages all distances to find rep. distance @ each node.
        3) Assigns a tau0 value based on depth and rep. distance.
        """
        if UnstructuredMesh.z is None:
            raise Exception('UnstructuredMesh must have z values.')
        if UnstructuredMesh.SpatialReference is None:
            raise Exception('UnstructuredMesh must have a spatial reference.')
        x = UnstructuredMesh.get_x(SpatialReference=3395)
        y = UnstructuredMesh.get_y(SpatialReference=3395)
        tri = Triangulation(x, y, UnstructuredMesh.elements)
        neighbors = defaultdict(set)
        for simplex in tri.triangles:
            for i, j in permutations(simplex, 2):
                neighbors[i].add(j)
        points = [tuple(p) for p in np.vstack([x, y]).T]
        values = np.full(UnstructuredMesh.z.shape, default_value)
        for k, v in neighbors.items():
            x0, y0 = points[k]
            distances = list()
            for idx in v:
                x1, y1 = points[idx]
                distances.append(np.sqrt((x0 - x1)**2 + (y0 - y1)**2))
            distance = np.mean(distances)
            if distance >= threshold_distance:
                if UnstructuredMesh.z[k] >= threshold_depth:
                    values[k] = shallow_tau0
                else:
                    values[k] = deep_tau0
        return values, default_value

    @property
    def TAU0(self):
        return self._TAU0

    @property
    def _TAU0(self):
        return self.__TAU0

    @_TAU0.setter
    def _TAU0(self, TAU0):
        self.__TAU0 = TAU0
