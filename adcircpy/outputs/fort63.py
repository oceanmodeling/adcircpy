import numpy as np
import pathlib
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
# from matplotlib.widgets import Slider
from copy import deepcopy
from adcircpy.mesh import UnstructuredMesh


class ElevationSurface(UnstructuredMesh):

    def __init__(
        self,
        vertices,
        elements,
        zeta,
        crs=None,
    ):
        super().__init__(vertices, elements, crs)
        self._zeta = zeta

    def write_2dm(self, path, index):
        path = pathlib.Path(path)
        triangles = self.triangulation.triangles
        f = "MESH2D\n"
        for i in range(triangles.shape[0]):
            f += f"E3T {i + 1} "
            for j in range(len(triangles[i, :])):
                f += f"{triangles[i, j]+1} "
            f += "\n"
        for i in range(self.xy.shape[0]):
            f += f"ND {i + 1} "
            f += f"{self.x[i]:<.16E} "
            f += f"{self.y[i]:<.16E} "
            f += f"{self.values[index, i]:<.16E}\n"
        with open(path.absolute(), 'w') as h:
            h.write(f)

    def make_movie(self, start_index=0, end_index=-1):

        fig = plt.figure(figsize=(18.5, 10))
        # mng = plt.get_current_fig_manager()
        # mng.resize(*mng.window.maxsize())
        axes = fig.add_subplot()
        fig.set_tight_layout(True)

        tri = self.triangulation

        mask = self.values.mask[0]
        tri.set_mask(np.any(mask[tri.triangles], axis=1))
        neg = np.mean(self.values[np.where(self.values[start_index:end_index] < 0)])
        pos = np.mean(self.values[np.where(self.values[start_index:end_index] > 0)])
        vmin = np.min([neg, -pos])
        vmax = np.max([-neg, pos])
        ax = axes.tripcolor(
            tri,
            self.values[0, :].data,
            # levels=256,
            cmap='RdBu_r',
            # norm=norm,
            vmin=vmin,
            # vmin=-0.3,
            vmax=vmax,
            # vmax=0.3,
            # **kwargs
            shading='gouraud'
            )
        plt.title(self.nc['time'][0])
        plt.gca().axis('scaled')
        plt.colorbar(ax)

        def update(i):
            mask = self.values.mask[i]
            tri.set_mask(np.any(mask[tri.triangles], axis=1))
            plt.title(f'NSPOOL {i}')
            ax.set_array(self.values[i, :].data)
            return ax, axes

        anim = FuncAnimation(
            fig,
            update,
            frames=np.arange(start_index, end_index),
            interval=200
            )
        plt.show()
        # anim.save(
        #     '/home/jreniel/ADCIRC_Maria2017_GAHM_sigma_null.gif',
        #     writer='imagemagick')
        return anim

    def make_plot(
        self,
        index,
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
        values = self.values[index, :]
        if axes is None:
            axes = plt.figure(figsize=figsize).add_subplot(111)
        if vmin is None:
            vmin = np.min(values)
        if vmax is None:
            vmax = np.max(values)
        tri = deepcopy(self.triangulation)
        tri_mask = np.any(values.mask[tri.triangles], axis=1)
        tri.set_mask(tri_mask)
        ax = axes.tricontourf(
            tri,
            values,
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
    def open(cls, path, crs=None, fort14=None):
        try:
            cls = cls.from_netcdf(path, crs)
            cls.nc = Dataset(path)
            return cls
        except OSError:
            return cls.from_ascii(path, fort14, crs)

    @classmethod
    def from_netcdf(cls, path, crs=None):
        nc = Dataset(path)
        cls._certify_netcdf_fort63_file(nc)
        zeta = np.ma.masked_equal(
            nc['zeta'][:], nc['zeta']._FillValue)
        return cls(np.vstack([nc['x'][:], nc['y'][:]]).T,
                   nc['element'][:]-1,
                   zeta,
                   crs)

    @classmethod
    def from_ascii(cls, path, fort14, crs=None):
        raise NotImplementedError
        # fort14 = cls.parse_gr3(fort14)
        # with open(path, 'r') as f:
        #     line = f.readline()
        #     line = f.readline().split()
        #     NP = int(line[1])
        #     line = f.readline()
        #     values = list()
        #     for i in range(NP):
        #         values.append(float(f.readline().split()[1]))
        #     for i in range(NP):
        #         try:
        #             .append(float(f.readline().split()[1]))
        #         except IndexError:
        #             .append(-99999.)
        # values = np.ma.masked_equal(values, -99999.)
        #  = np.ma.masked_equal(, -99999.)
        # return cls(
        #     np.vstack([fort14.pop('x'), fort14.pop('y')]).T,
        #     fort14.pop('elements'),
        #     values,
        #     ,
        #     crs)

    @staticmethod
    def _certify_netcdf_fort63_file(nc):
        if ('zeta' not in nc.variables.keys()
                and 'adcirc_mesh' not in nc.variables.keys()):
            raise Exception('Not a fort63 file!')

    @property
    def values(self):
        return self.zeta

    @property
    def zeta(self):
        return self._zeta

    @property
    def _zeta(self):
        return self.__zeta

    @_zeta.setter
    def _zeta(self, zeta):
        self.__zeta = zeta


# alias
Fort63 = ElevationSurface
