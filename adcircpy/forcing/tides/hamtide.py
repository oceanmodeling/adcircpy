from functools import lru_cache
from os import PathLike
from pathlib import Path

from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import griddata

from adcircpy.forcing.tides.dataset import TidalDataset


class HAMTIDE(TidalDataset):
    '''
    Taguchi, E., D. Stammer and W. Zahel (2010), Estimation of deep ocean tidal energy dissipation based on the high-resolution data-assimilative HAMTIDE model (to be submitted to J. Geophys. Res.).
    https://icdc.cen.uni-hamburg.de/en/hamtide.html
    '''

    CONSTITUENTS = ['S2', 'Q1', 'P1', 'O1', 'N2', 'M2', 'K2', 'K1']
    OPENDAP_URL = 'https://icdc.cen.uni-hamburg.de/thredds/dodsC/ftpthredds/hamtide/'

    def __init__(self, hamtide_dataset_directory: PathLike = None):
        if hamtide_dataset_directory is not None:
            hamtide_dataset_directory = Path(hamtide_dataset_directory)
            if len(list(hamtide_dataset_directory.glob('*.nc'))) == 0:
                raise FileNotFoundError(
                        f'no NetCDF files found at "{hamtide_dataset_directory}"')
        else:
            hamtide_dataset_directory = self.OPENDAP_URL

        super().__init__(hamtide_dataset_directory)

        datasets = {'elevation': {}, 'velocity': {}}
        for variable in datasets.keys():
            datasets[variable].update(
                    {
                        constituent: {'path': None, 'dataset': None}
                        for constituent in self.CONSTITUENTS
                    }
            )

        self.datasets = datasets

    def get_amplitude(
            self,
            constituent: str,
            vertices: np.ndarray
    ) -> np.ndarray:
        if not isinstance(vertices, np.ndarray):
            vertices = np.asarray(vertices)
        self._assert_vertices(vertices)
        return self._get_interpolation('elevation', 'AMPL', constituent,
                                       vertices)

    def get_phase(
            self,
            constituent: str,
            vertices: np.ndarray
    ) -> np.ndarray:
        if not isinstance(vertices, np.ndarray):
            vertices = np.asarray(vertices)
        self._assert_vertices(vertices)
        return self._get_interpolation('elevation', 'PHAS', constituent,
                                       vertices)

    @property
    @lru_cache(maxsize=1)
    def x(self) -> np.ndarray:
        return Dataset(self._prepend_path('k2.hamtide11a.nc'))['LON'][:].data

    @property
    @lru_cache(maxsize=1)
    def y(self) -> np.ndarray:
        return Dataset(self._prepend_path('k2.hamtide11a.nc'))['LAT'][:].data

    def _get_dataset(self, variable: str, constituent: str) -> Dataset:
        data = self.datasets[variable][constituent]

        dataset = data['dataset']
        if dataset is None:
            path = data['path']
            if path is None:
                if variable == 'elevation':
                    filename = f'{constituent.lower()}.hamtide11a.nc'
                elif variable == 'velocity':
                    filename = f'HAMcurrent11a_{constituent.lower()}.nc'
                else:
                    raise NotImplementedError(f'tidal variable "{variable}" '
                                              f'not implemented')

                path = self._prepend_path(filename)

            try:
                dataset = Dataset(path)
                if data['path'] is None:
                    self.datasets[variable][constituent]['path'] = path
                self.datasets[variable][constituent]['dataset'] = dataset
            except FileNotFoundError:
                raise FileNotFoundError(f'no dataset found for "{variable}" '
                                        f'"{constituent}" at "{path}"')

        return dataset

    def _get_interpolation(
            self,
            variable: str,
            netcdf_variable: str,
            constituent: str,
            vertices: np.ndarray
    ) -> np.ndarray:
        self._assert_vertices(vertices)

        xq = np.asarray([x + 360. if x < 0. else x
                         for x in vertices[:, 0]]).flatten()
        yq = vertices[:, 1].flatten()

        xidx = np.logical_and(self.x >= np.min(xq), self.x <= np.max(xq))
        yidx = np.logical_and(self.y >= np.min(yq), self.y <= np.max(yq))

        dataset = self._get_dataset(variable, constituent)

        xi, yi = np.meshgrid(self.x[xidx], self.y[yidx])
        xi = xi.flatten()
        yi = yi.flatten()
        zi = dataset[netcdf_variable][yidx, xidx].flatten()

        return griddata(
                (xi[~zi.mask], yi[~zi.mask]),
                zi[~zi.mask],
                (xq, yq),
                method='nearest'
        ) * 0.01

    def _prepend_path(self, filename: str) -> str:
        if self.path is None:
            path = self.path
        elif isinstance(self.path, Path):
            path = self.path / filename
        else:
            path = f'{self.path}/{filename}'
        return path
