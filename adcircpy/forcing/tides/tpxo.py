import os
from os import PathLike
from pathlib import Path

import appdirs
from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import griddata

from adcircpy.forcing.tides.dataset import TidalDataset

TPXO_ENVIRONMENT_VARIABLE = 'TPXO_NCFILE'
TPXO_FILENAME = 'h_tpxo9.v1.nc'


class TPXO(TidalDataset):
    DEFAULT_PATH = Path(appdirs.user_data_dir('tpxo')) / TPXO_FILENAME

    def __init__(self, tpxo_dataset_filename: PathLike = None):
        if tpxo_dataset_filename is None:
            tpxo_environment_variable = os.getenv(TPXO_ENVIRONMENT_VARIABLE)
            if tpxo_environment_variable is not None:
                tpxo_dataset_filename = tpxo_environment_variable
            else:
                tpxo_dataset_filename = self.DEFAULT_PATH

        super().__init__(tpxo_dataset_filename)

        if self.path is not None:
            self.dataset = Dataset(self.path)
        else:
            raise FileNotFoundError(
                '\n'.join(
                    [
                        f'No TPXO file found at "{self.path}".',
                        'New users will need to register and request a copy of '
                        f'the TPXO9 NetCDF file (specifically `{TPXO_FILENAME}`) '
                        'from the authors at https://www.tpxo.net.',
                        'Once you obtain `h_tpxo9.v1.nc`, you can follow one of the following options: ',
                        f'1) copy or symlink the file to "{self.path}"',
                        f'2) set the environment variable `{TPXO_ENVIRONMENT_VARIABLE}` to point to the file',
                    ]
                )
            )

    def get_amplitude(self, constituent: str, vertices: np.ndarray) -> np.ndarray:
        if not isinstance(vertices, np.ndarray):
            vertices = np.asarray(vertices)
        self._assert_vertices(vertices)
        return self._get_interpolation(self.ha, constituent, vertices)

    def get_phase(self, constituent: str, vertices: np.ndarray) -> np.ndarray:
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

    @property
    def constituents(self):
        if not hasattr(self, '_constituents'):
            self._constituents = [
                c.capitalize()
                for c in self.dataset['con'][:]
                .astype('|S1')
                .tostring()
                .decode('utf-8')
                .split()
            ]
        return self._constituents

    def _get_interpolation(
        self, tpxo_array: np.ndarray, constituent: str, vertices: np.ndarray
    ):
        """
        `tpxo_index_key` is either `ha` or `hp` based on the keys used
        internally in the TPXO NetCDF file.
        """

        self._assert_vertices(vertices)
        constituents = list(map(lambda x: x.lower(), self.constituents))
        constituent = constituents.index(constituent.lower())
        array = tpxo_array[constituent, :, :].flatten()
        _x = np.asarray([x + 360.0 for x in vertices[:, 0] if x < 0]).flatten()
        _y = vertices[:, 1].flatten()
        x, y = np.meshgrid(self.x, self.y, indexing='ij')
        x = x.flatten()
        y = y.flatten()
        dx = np.mean(np.diff(self.x))
        dy = np.mean(np.diff(self.y))

        # buffer the bbox by 2 difference units
        _idx = np.where(
            np.logical_and(
                np.logical_and(x >= np.min(_x) - 2 * dx, x <= np.max(_x) + 2 * dx),
                np.logical_and(y >= np.min(_y) - 2 * dy, y <= np.max(_y) + 2 * dy),
            )
        )

        # "method" can be 'spline' or any string accepted by griddata()'s method kwarg.
        values = griddata(
            (x[_idx], y[_idx]), array[_idx], (_x, _y), method='linear', fill_value=np.nan,
        )
        nan_idxs = np.where(np.isnan(values))
        values[nan_idxs] = griddata(
            (x[_idx], y[_idx]), array[_idx], (_x[nan_idxs], _y[nan_idxs]), method='nearest',
        )
        return values
