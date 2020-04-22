from scipy.interpolate import griddata
import numpy as np
import os
import sys
from netCDF4 import Dataset
import pathlib


class TPXO:
    """
    Egbert, Gary D., and Svetlana Y. Erofeeva. "Efficient inverse modeling of barotropic ocean tides." Journal of Atmospheric and Oceanic Technology 19.2 (2002): 183-204
    """

    def __init__(self):
        file = os.getenv('TPXO_NCFILE')
        if file is not None:
            self._nc = Dataset(file)
            return

        else:
            prefix = "/".join(sys.executable.split('/')[:-2])
            file = pathlib.Path(prefix) / 'lib/h_tpxo9.v1.nc'
        if isinstance(file, pathlib.Path):
            if file.is_file():
                self._nc = Dataset(file)
                return

        else:
            msg = "No TPXO file found. You need to register and request a "
            msg += "copy of the TPXO9 netcdf file (specifically h_tpxo9.v1.nc)"
            msg += " from the authors at https://www.tpxo.net. Once you obtain"
            msg += " this copy, set the environment variable TPXO_NCFILE "
            msg += "to point to the path of the h_tpxo9.v1.nc file."
            raise FileNotFoundError(msg)

    def __call__(self, constituent, vertices):
        """
        "method" can be 'spline' or any string accepted by griddata()'s method
        kwarg.
        """
        amp = self.get_amplitude(constituent, vertices)
        phase = self.get_phase(constituent, vertices)
        return (amp, phase)

    def get_amplitude(self, constituent, vertices):
        """
        "method" can be 'spline' or any string accepted by griddata()'s method
        kwarg.
        """
        vertices = np.asarray(vertices)
        self._assert_vertices(vertices)
        return self._get_interpolation(self.ha, constituent, vertices)

    def get_phase(self, constituent, vertices):
        vertices = np.asarray(vertices)
        self._assert_vertices(vertices)
        return self._get_interpolation(self.hp, constituent, vertices)

    @property
    def x(self):
        return self.nc['lon_z'][:, 0].data

    @property
    def y(self):
        return self.nc['lat_z'][0, :].data

    @property
    def ha(self):
        return self.nc['ha'][:]

    @property
    def hp(self):
        return self.nc['hp'][:]

    @property
    def nc(self):
        return self._nc

    @property
    def tpxo_constituents(self):
        return ['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1', 'Mm', 'Mf',
                'M4', 'MN4', 'MS4', '2N2', 'S1']

    def _get_interpolation(self, tpxo_array, constituent, vertices):
        """
        tpxo_index_key is either 'ha' or 'hp' based on the keys used internally
        on the TPXO netcdf file.
        """
        constituent = self.tpxo_constituents.index(constituent)
        array = tpxo_array[constituent, :, :].flatten()
        _x = np.asarray([x + 360. for x in vertices[:, 0] if x < 0]).flatten()
        _y = vertices[:, 1].flatten()
        x, y = np.meshgrid(self.x, self.y, indexing='ij')
        x = x.flatten()
        y = y.flatten()
        dx = np.mean(np.diff(self.x))
        dy = np.mean(np.diff(self.y))
        _idx = np.where(
            np.logical_and(  # buffer the bbox by 2 difference units
                np.logical_and(x >= np.min(_x)-2*dx, x <= np.max(_x)+2*dx),
                np.logical_and(y >= np.min(_y)-2*dy, y <= np.max(_y)+2*dy)))
        return griddata(
            (x[_idx], y[_idx]), array[_idx], (_x, _y), method='nearest')

    def _assert_vertices(self, vertices):
        msg = "vertices must be of shape M x 2"
        assert vertices.shape[1] == 2, msg

    @property
    def _nc(self):
        return self.__nc

    @_nc.setter
    def _nc(self, nc):
        self.__nc = nc


def install():
    import shutil
    prefix = "/".join(sys.executable.split('/')[:-2])
    shutil.copyfile(sys.argv[1], pathlib.Path(prefix) / 'lib/h_tpxo9.v1.nc')
