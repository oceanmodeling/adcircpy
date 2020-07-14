import numpy as np
import pathlib
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from adcircpy.mesh import grd
from adcircpy.mesh import EuclideanMesh2D


class Maxele(EuclideanMesh2D):

    def __init__(
        self,
        nodes,
        elements,
        time_of_zeta_max=None,
        description='maxele',
        crs=None,
    ):
        super().__init__(**grd.euclidean_mesh({
            'nodes': nodes,
            'elements': elements,
            'description': description,
            'crs': crs
            }))
        self._zeta_max = np.ma.masked_equal(super().values, 99999.0)
        self._time_of_zeta_max = time_of_zeta_max

    def make_plot(
        self,
        axes=None,
        vmin=None,
        vmax=None,
        cmap='jet',
        levels=256,
        show=False,
        title=None,
        figsize=None,
        extent=None,
        cbar_label=None,
        **kwargs
    ):
        if axes is None:
            axes = plt.figure(figsize=figsize).add_subplot(111)
        if vmin is None:
            vmin = np.min(self.values)
        if vmax is None:
            vmax = np.max(self.values)
        tri = self.triangulation
        if np.any(self.values.mask):
            tri_mask = np.any(self.values.mask[tri.triangles], axis=1)
            tri.set_mask(tri_mask)
        ax = axes.tricontourf(
            tri,
            self.values.data,
            levels=levels,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            **kwargs
            )
        plt.colorbar(ax, cmap=cmap)
        axes.axis('scaled')
        if extent is not None:
            axes.axis(extent)
        if show is True:
            plt.show()
        return axes

    @classmethod
    def open(cls, path, crs=None):
        path = str(pathlib.Path(path).absolute())
        try:
            return cls.from_netcdf(path, crs)
        except OSError:
            return cls.from_ascii(path, crs)

    @classmethod
    def from_netcdf(cls, path, crs=None):
        nc = Dataset(path)
        cls._certify_netcdf_maxele_file(nc)
        zeta_max = np.ma.masked_equal(
            nc['zeta_max'][:], nc['zeta_max']._FillValue)
        time_of_zeta_max = np.ma.masked_equal(
            nc['time_of_zeta_max'][:], nc['time_of_zeta_max']._FillValue)
        xyz = np.vstack([nc['x'][:], nc['y'][:], zeta_max]).T
        nodes = {i: ((x, y), z) for i, (x, y, z) in enumerate(xyz)}
        elements = {
            i: (e0-1, e1-1, e2-1) for i, (e0, e1, e2) in enumerate(nc['element'][:])
        }
        return cls(nodes, elements,
                   # zeta_max,
                   time_of_zeta_max,
                   'maxele',
                   crs
                   )

    @classmethod
    def from_ascii(cls, path, fort14, crs=None):
        fort14 = cls.parse_gr3(fort14)
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
        return cls(
            np.vstack([fort14.pop('x'), fort14.pop('y')]).T,
            fort14.pop('elements'),
            values,
            time_of_zeta_max,
            crs)

    @staticmethod
    def _certify_netcdf_maxele_file(nc):
        if ('zeta_max' not in nc.variables.keys()
                and 'time_of_zeta_max' not in nc.variables.keys()):
            raise Exception('Not a maxele file!')

    @property
    def values(self):
        return self.zeta_max

    @property
    def zeta_max(self):
        return self._zeta_max

    @property
    def time_of_zeta_max(self):
        return self._time_of_zeta_max

    @property
    def _zeta_max(self):
        return self.__zeta_max

    @property
    def _time_of_zeta_max(self):
        return self.__time_of_zeta_max

    @_zeta_max.setter
    def _zeta_max(self, zeta_max):
        self.__zeta_max = zeta_max

    @_time_of_zeta_max.setter
    def _time_of_zeta_max(self, time_of_zeta_max):
        self.add_attribute('time_of_zeta_max')
        self.set_attribute('time_of_zeta_max', time_of_zeta_max)
