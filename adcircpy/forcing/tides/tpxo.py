import os
from os import PathLike
from pathlib import Path
import sys

from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import griddata

from adcircpy.forcing.tides.dataset import TidalDataset


class TPXO(TidalDataset):
    CONSTITUENTS = ['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1', 'Mm', 'Mf',
                    'M4', 'MN4', 'MS4', '2N2', 'S1']
    DEFAULT_PATH = Path(sys.executable).parent.parent / 'lib' / 'h_tpxo9.v1.nc'

    def __init__(self, path: PathLike = None):
        if path is None:
            environment_variable = os.getenv('TPXO_NCFILE')
            if environment_variable is not None:
                path = environment_variable
            else:
                path = self.DEFAULT_PATH

        super().__init__(path)

        if self.path is not None:
            self.dataset = Dataset(self.path)
        else:
            raise FileNotFoundError(
                    'No TPXO file found.\n'
                    'New users will need to register and request a copy of '
                    'the TPXO9 NetCDF file (specifically `h_tpxo9.v1.nc`) '
                    'from the authors at https://www.tpxo.net.\n'
                    'Once you obtain this copy, you can set '
                    'the environment variable `TPXO_NCFILE` to point to '
                    'the path of the `h_tpxo9.v1.nc` file. '
                    f'Alternatively, you can symlink this file manually to '
                    f'"{self.DEFAULT_PATH}"'
            )

    def __call__(self, constituent: str, vertices) -> (np.ndarray, np.ndarray):
        """
        "method" can be 'spline' or any string accepted by griddata()'s method
        kwarg.
        """
        return self.get_amplitude(constituent, vertices), \
               self.get_phase(constituent, vertices)

    def get_amplitude(
            self,
            constituent: str,
            vertices: np.ndarray
    ) -> np.ndarray:
        if not isinstance(vertices, np.ndarray):
            vertices = np.asarray(vertices)
        self._assert_vertices(vertices)
        return self._get_interpolation(self.ha, constituent, vertices)

    def get_phase(
            self,
            constituent: str,
            vertices: np.ndarray
    ) -> np.ndarray:
        if not isinstance(vertices, np.ndarray):
            vertices = np.asarray(vertices)
        self._assert_vertices(vertices)
        return self._get_interpolation(self.hp, constituent, vertices)

    @property
    def x(self) -> np.ndarray:
        return self.dataset['lon_z'][:, 0].data

    @property
    def y(self) -> np.ndarray:
        return self.dataset['lat_z'][0, :].data

    @property
    def ha(self) -> np.ndarray:
        return self.dataset['ha'][:]

    @property
    def hp(self) -> np.ndarray:
        return self.dataset['hp'][:]

    def _get_interpolation(self, tpxo_array: np.ndarray, constituent: str,
                           vertices: np.ndarray):
        """
        `tpxo_index_key` is either `ha` or `hp` based on the keys used
        internally in the TPXO NetCDF file.
        """

        self._assert_vertices(vertices)

        constituent = self.CONSTITUENTS.index(constituent)
        array = tpxo_array[constituent, :, :].flatten()
        _x = np.asarray([x + 360. for x in vertices[:, 0] if x < 0]).flatten()
        _y = vertices[:, 1].flatten()
        x, y = np.meshgrid(self.x, self.y, indexing='ij')
        x = x.flatten()
        y = y.flatten()
        dx = np.mean(np.diff(self.x))
        dy = np.mean(np.diff(self.y))

        # buffer the bbox by 2 difference units
        _idx = np.where(
                np.logical_and(
                        np.logical_and(
                                x >= np.min(_x) - 2 * dx,
                                x <= np.max(_x) + 2 * dx
                        ),
                        np.logical_and(
                                y >= np.min(_y) - 2 * dy,
                                y <= np.max(_y) + 2 * dy
                        )
                )
        )

        # "method" can be 'spline' or any string accepted by griddata()'s method kwarg.
        return griddata(
                (x[_idx], y[_idx]),
                array[_idx],
                (_x, _y),
                method='nearest'
        )
