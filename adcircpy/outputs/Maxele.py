import warnings
import numpy as np
from netCDF4 import Dataset
from AdcircPy.Model import AdcircMesh
from AdcircPy.Outputs.ScalarSurfaceExtrema import ScalarSurfaceExtrema


class Maxele(ScalarSurfaceExtrema):
    def __init__(self, xy, values, elements, SpatialReference, vertical_datum,
                 model_times=None, node_id=None, element_id=None,
                 **Boundaries):
        super(Maxele, self).__init__(
            xy, values, elements, SpatialReference, vertical_datum,
            model_times, node_id, element_id, **Boundaries)

    def convert_datum(self, datum_mesh, target_datum, operator='plus'):
        """
        The way vertical datums are implemented is not self consistent.
        Vertical datum grids have to be handled with care until
        some standards can be set. Meshes should actually be referenced
        to an equipotential surface, and not to space variable datums.
        This is because at time t0 in the model run, the water surface
        elevation must sit on an equipotential surface. Datums derived from
        tidal measurements are known to not be at the same geopotentials.
        """
        assert operator in ['plus', 'minus']
        fort14 = AdcircMesh.parse_fort14(datum_mesh)
        x = np.asarray(fort14['x'])
        y = np.asarray(fort14['y'])
        xy = np.vstack([x, y]).T
        try:
            np.testing.assert_array_equal(xy, self.xy)
        except AssertionError:
            assert xy.shape == self.xy.shape
            warnings.warn(
                "Coordinates of datum mesh are not identical to mesh, but "
                + "sizes do match. Continuing...")

        values = np.asarray(fort14['z'])
        if operator == 'plus':
            self._z = self.z + values
        else:
            self._z = self.z - values
        self._vertical_datum = target_datum

    @classmethod
    def from_netcdf(cls, path, fort14=None, vertical_datum='LMSL',
                    SpatialReference=4326):
        nc = Dataset(path)
        if 'zeta_max' not in nc.variables.keys():
            raise Exception('Not a maxele file!')
        fort14 = cls.__init_fort14(fort14, SpatialReference, vertical_datum)
        xy = np.vstack([nc['x'][:], nc['y'][:]]).T
        return cls(xy,
                   nc['zeta_max'][:],
                   nc['element'][:]-1,
                   SpatialReference,
                   vertical_datum,
                   nc['time_of_zeta_max'])

    @classmethod
    def from_ascii(cls, path, fort14, vertical_datum='LMSL',
                   SpatialReference=4326, datum_grid=None):
        fort14 = cls.__init_fort14(fort14)
        with open(path, 'r') as f:
            line = f.readline()
            line = f.readline().split()
            NP = int(line[1])
            line = f.readline()
            values = list()
            for i in range(NP):
                values.append(float(f.readline().split()[1]))
            time_of_zeta_max = list()
            for i in range(NP):
                try:
                    time_of_zeta_max.append(float(f.readline().split()[1]))
                except IndexError:
                    time_of_zeta_max.append(-99999.)

        values = np.ma.masked_equal(values, -99999.)
        time_of_zeta_max = np.ma.masked_equal(time_of_zeta_max, -99999.)
        return cls(fort14.x,
                   fort14.y,
                   fort14.elements,
                   values,
                   time_of_zeta_max,
                   SpatialReference=SpatialReference,
                   vertical_datum=vertical_datum,
                   ocean_boundaries=fort14.ocean_boundaries,
                   land_boundaries=fort14.land_boundaries,
                   inner_boundaries=fort14.inner_boundaries,
                   weir_boundaries=fort14.weir_boundaries,
                   inflow_boundaries=fort14.inflow_boundaries,
                   outflow_boundaries=fort14.outflow_boundaries,
                   culvert_boundaries=fort14.culvert_boundaries,
                   datum_grid=datum_grid)

    @staticmethod
    def __init_fort14(fort14, SpatialReference, vertical_datum):
        if fort14 is not None:
            if isinstance(fort14, AdcircMesh) is False:
                return AdcircMesh(fort14, SpatialReference, vertical_datum)
