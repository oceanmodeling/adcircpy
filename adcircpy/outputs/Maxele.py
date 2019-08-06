import numpy as np
from pathlib import Path
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from adcircpy.mesh import UnstructuredMesh


class Maxele(UnstructuredMesh):

    def __init__(
        self,
        vertices,
        elements,
        zeta_max,
        time_of_zeta_max,
        SpatialReference=None,
    ):
        super(Maxele, self).__init__(
            vertices, elements, zeta_max, SpatialReference)
        self.__time_of_zeta_max = time_of_zeta_max

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
        ax = axes.tricontourf(
            self.mpl_tri, self.values, levels=levels,
            cmap=cmap, vmin=vmin, vmax=vmax, **kwargs)
        plt.colorbar(ax, cmap=cmap)
        axes.axis('scaled')
        if extent is not None:
            axes.axis(extent)
        if show is True:
            plt.show()
        return axes

    @classmethod
    def from_netcdf(cls, path, SpatialReference=None):
        nc = Dataset(str(Path(path)))
        if 'zeta_max' not in nc.variables.keys():
            raise Exception('Not a maxele file!')
        values = np.ma.masked_equal(
            nc['zeta_max'][:], nc['zeta_max']._FillValue)
        return cls(np.vstack([nc['x'][:], nc['y'][:]]).T,
                   nc['element'][:]-1,
                   values,
                   SpatialReference)

    @classmethod
    def from_ascii(cls, path, fort14, SpatialReference=None):
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
            fort14.pop('elements'), values, time_of_zeta_max, SpatialReference)

    @property
    def time_of_zeta_max(self):
        return self.__time_of_zeta_max
