import abc
from functools import lru_cache
import pathlib

from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from netCDF4 import Dataset
import numpy as np
from pyproj import CRS

from adcircpy.figures import figure
from adcircpy.mesh.parsers import sms2dm

# class OutputVariable(Enum):
#     FORT63 = 'zeta'


class SurfaceOutput(metaclass=abc.ABCMeta):
    # change this for __types__

    _physical_variables = {
        'fort.63': 'zeta',
        'maxele': 'zeta_max',
        'time_of_maxele': 'time_of_zeta_max',
    }

    def __init__(self, path, crs=None):
        self._path = path
        self._crs = crs

    def export(self, path, overwrite=False):
        coords = {i + 1: (self.x[i], self.y[i]) for i in range(len(self.values))}
        values = self.values.filled(99999.0)
        nodes = {id: ((x, y), values[i]) for i, (id, (x, y)) in enumerate(coords.items())}
        triangles = {id + 1: tuple(e) for id, e in enumerate(self.triangles + 1)}
        sms2dm.writer({'ND': nodes, 'E3T': triangles}, path, overwrite)

    @figure
    def triplot(self, *args, axes=None, color='k', linewidth=0.1, **kwargs):
        plt.triplot(self.x, self.y, self.triangles, color=color, linewidth=linewidth)
        return axes

    @figure
    def tricontourf(self, *args, axes=None, **kwargs):
        if np.any(self.values.mask):
            self.triangulation.set_mask(
                np.any(self.values.mask[self.triangulation.triangles], axis=1)
            )
        _ax = plt.tricontourf(
            self.triangulation,
            self.values,
            cmap=kwargs.get('cmap', self._cmap),
            levels=kwargs.get('levels', self._levels),
            vmin=kwargs.get('vmin', np.min(self.values)),
            vmax=kwargs.get('vmax', np.max(self.values)),
        )
        self.triangulation.set_mask(None)
        if kwargs.get('cbar') is not None:
            plt.colorbar(_ax)
        plt.gca().axis('scaled')
        return axes

    @property
    @lru_cache(maxsize=None)
    def x(self):
        if isinstance(self._ptr, Dataset):
            return self._ptr['x'][:].data
        else:
            raise NotImplementedError('ascii')

    @property
    @lru_cache(maxsize=None)
    def y(self):
        if isinstance(self._ptr, Dataset):
            return self._ptr['y'][:].data
        else:
            raise NotImplementedError('ascii')

    @property
    @lru_cache(maxsize=None)
    def triangles(self):
        if isinstance(self._ptr, Dataset):
            return self._ptr['element'][:].data - 1
        else:
            raise NotImplementedError('ascii')

    @property
    @lru_cache(maxsize=None)
    def triangulation(self):
        return Triangulation(self.x, self.y, triangles=self.triangles)

    @property
    def values(self):
        return self._values

    @property
    def crs(self):
        return self.__crs

    @property
    def _path(self):
        return self.__path

    @_path.setter
    def _path(self, path):
        path = pathlib.Path(path)

        if not path.is_file():
            raise IOError(f'File not found: {str(path)}')
        self.__path = path

    @property
    def _crs(self):
        return self.__crs

    @_crs.setter
    def _crs(self, crs):
        if crs is not None:
            crs = CRS.from_user_input(crs)
        self.__crs = crs

    @property
    def _cmap(self):
        return None

    @property
    def _levels(self):
        return None

    @property
    @lru_cache(maxsize=None)
    def _physical_variable(self):
        return self._physical_variables[self._filetype]

    @property
    @abc.abstractmethod
    def _filetype(self):
        """Subclass must implement _filetype attribute."""

    @property
    def _values(self):
        try:
            return self.__values
        except AttributeError:
            return self._ptr[self._physical_variable][:]

    @_values.setter
    def _values(self, values):
        self.__values = values

    @property
    @lru_cache(maxsize=None)
    def _ptr(self):
        try:
            nc = Dataset(self._path)
            msg = 'NetCDF file provided is not a surface output.'
            assert 'adcirc_mesh' in nc.variables, msg
            msg = f'"{self._physical_variable}" variable not found in file: '
            msg += f'{self._path}, '
            msg += f'therefore, this is not a {self._filetype} file.'
            assert self._physical_variable in nc.variables, msg
            return nc
        except OSError as e:
            if e.errno == -51:
                if self._is_ascii:
                    raise NotImplementedError('ASCII outputs are not implemented.')
                    values = self._get_ascii_values(0)
                    return {f'{self._physical_variable}': values}
            raise e

    @property
    def _is_ascii(self):
        try:
            with open(self._path, 'r') as f:
                f.readline()
                if 'FileFmtVersion' in f.readline():
                    return True
            return False
        except:  # noqa: E722
            return False

    @property
    def _is_netcdf(self):
        try:
            Dataset(self._path)
        except:  # noqa: E722
            return False


class SurfaceOutputTimeseries(SurfaceOutput):
    def __init__(self, path, crs=None, index=0):
        super().__init__(path, crs)
        self.index = index

    def __iter__(self):
        for i in range(len(self)):
            yield self._ptr[self._physical_variable][i, :]

    def __next__(self):
        if self.index < len(self):
            self.index += 1
            return self.values
        else:
            raise StopIteration

    def __len__(self):
        return self._ptr.dimensions['time'].size

    @property
    def index(self):
        return self.__index

    @index.setter
    def index(self, index):
        if isinstance(self._ptr, Dataset):
            self._values = self._ptr[self._physical_variable][index, :]
        else:
            raise NotImplementedError('ascii')
        self.__index = index

    @property
    @abc.abstractmethod
    def animation(self):
        """Subclass must implement animation method."""


class ScalarSurfaceOutputTimeseries(SurfaceOutputTimeseries):
    def animation(self, save=False, fps=3, start_frame=0, end_frame=-1, **kwargs):

        fig = plt.figure(figsize=kwargs.get('figsize'))
        ax = fig.add_subplot(111)
        plt.tight_layout(pad=2)
        _oi = self.index

        def animate(i):
            self.index = i
            _ax = fig.get_axes()
            ax.clear()
            if len(_ax) > 1:
                cax = _ax[1]
                cax.cla()
            else:
                cax = None
            if np.any(self.values.mask):
                self.triangulation.set_mask(
                    np.any(self.values.mask[self.triangulation.triangles], axis=1)
                )

            if kwargs.get('elements', False):
                ax.triplot(self.triangulation, color='k', linewidth=0.7)

            _ax = ax.tricontourf(
                self.triangulation,
                self.values,
                cmap=kwargs.get('cmap', self._cmap),
                levels=kwargs.get('levels', self._levels),
            )
            ax.set_ylim(ymin=kwargs.get('ymin'), ymax=kwargs.get('ymax'), auto=True)
            ax.set_xlim(xmin=kwargs.get('xmin'), xmax=kwargs.get('xmax'), auto=True)
            # ax.set_title(dates[i].strftime('%b %d, %Y %H:%M'))
            ax.set_xlabel('Longitude (°E)')
            ax.set_ylabel('Latitude (°N)')
            # cbar = fig.colorbar(_ax, cax=cax, format='%.1f')
            # cbar.ax.set_ylabel('UNITS', rotation=90)

        end_frame = end_frame % len(self) if end_frame < 0 else end_frame
        start_frame = start_frame % len(self) if start_frame < 0 else start_frame
        frames = range(start_frame, end_frame)
        anim = FuncAnimation(fig, animate, frames, blit=False)

        if save:
            anim.save(pathlib.Path(save), writer='imagemagick', fps=fps)

        self.index = _oi

        if kwargs.get('show', False):
            plt.show()

        return anim
