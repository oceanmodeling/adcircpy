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
        zi = tpxo_array[constituent, :, :].flatten()
        xo = np.asarray([x + 360.0 for x in vertices[:, 0] if x < 0]).flatten()
        yo = vertices[:, 1].flatten()
        xi, yi = np.meshgrid(self.x, self.y, indexing='ij')
        xi = xi.flatten()
        yi = yi.flatten()
        dx = np.mean(np.diff(self.x))
        dy = np.mean(np.diff(self.y))
        # buffer window by 2 pixel units
        mask1 = np.logical_and(
            np.logical_and(xi >= np.min(xo) - 2 * dx, xi <= np.max(xo) + 2 * dx),
            np.logical_and(yi >= np.min(yo) - 2 * dy, yi <= np.max(yo) + 2 * dy),
        )
        # remove junk values from input array
        mask2 = np.ma.masked_where(zi != 0.0, zi)
        iidx = np.where(np.logical_and(mask1, mask2))
        values = griddata(
            (xi[iidx], yi[iidx]), zi[iidx], (xo, yo), method='linear', fill_value=np.nan,
        )
        nan_idxs = np.where(np.isnan(values))
        values[nan_idxs] = griddata(
            (xi[iidx], yi[iidx]), zi[iidx], (xo[nan_idxs], yo[nan_idxs]), method='nearest',
        )
        return values
