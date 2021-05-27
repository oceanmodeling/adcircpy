from os import PathLike
from pathlib import Path

from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import griddata

from adcircpy.forcing.tides.dataset import TidalDataset


class HAMTIDE(TidalDataset):
    """
    Taguchi, E., Stammer, D., & Zahel, W. (2010). Estimation of deep ocean tidal energy dissipation based on the high-resolution data-assimilative HAMTIDE model. J. geophys. Res.
    https://icdc.cen.uni-hamburg.de/en/hamtide.html
    """

    OPENDAP_URL = 'https://icdc.cen.uni-hamburg.de/thredds/dodsC/ftpthredds/hamtide/'

    def __init__(self, hamtide_dataset_directory: PathLike = None):
        if hamtide_dataset_directory is None:
            hamtide_dataset_directory = self.OPENDAP_URL
        else:
            try:
                if Path(hamtide_dataset_directory).exists():
                    hamtide_dataset_directory = Path(hamtide_dataset_directory)
                    if len(list(hamtide_dataset_directory.glob('*.nc'))) == 0:
                        raise FileNotFoundError(
                            f'no NetCDF files found at ' f'"{hamtide_dataset_directory}"'
                        )
            except OSError:
                raise ValueError('given resource must be a local path')

        super().__init__(hamtide_dataset_directory)

        datasets = {'elevation': {}, 'velocity': {}}
        for variable in datasets.keys():
            datasets[variable].update(
                {
                    constituent.lower(): {'path': None, 'dataset': None}
                    for constituent in self.constituents
                }
            )

        self.datasets = datasets

    def get_amplitude(self, constituent: str, vertices: np.ndarray) -> np.ndarray:
        if not isinstance(vertices, np.ndarray):
            vertices = np.asarray(vertices)
        self._assert_vertices(vertices)
        return self._get_interpolation('elevation', 'AMPL', constituent, vertices) * 0.01

    def get_phase(self, constituent: str, vertices: np.ndarray) -> np.ndarray:
        if not isinstance(vertices, np.ndarray):
            vertices = np.asarray(vertices)
        self._assert_vertices(vertices)
        return self._get_interpolation('elevation', 'PHAS', constituent, vertices)

    @property
    def x(self) -> np.ndarray:
        if not hasattr(self, '_x'):
            self._x = Dataset(self._prepend_path('k2.hamtide11a.nc'))['LON'][:].data
        return self._x

    @property
    def y(self) -> np.ndarray:
        if not hasattr(self, '_y'):
            self._y = Dataset(self._prepend_path('k2.hamtide11a.nc'))['LAT'][:].data
        return self._y

    @property
    def constituents(self):
        return ['S2', 'Q1', 'P1', 'O1', 'N2', 'M2', 'K2', 'K1']

    def _get_dataset(self, variable: str, constituent: str) -> Dataset:
        data = self.datasets[variable][constituent.lower()]

        dataset = data['dataset']
        if dataset is None:
            path = data['path']
            if path is None:
                if variable == 'elevation':
                    filename = f'{constituent.lower()}.hamtide11a.nc'
                elif variable == 'velocity':
                    filename = f'HAMcurrent11a_{constituent.lower()}.nc'
                else:
                    raise NotImplementedError(
                        f'tidal variable "{variable}" ' f'not implemented'
                    )

                path = self._prepend_path(filename)

            try:
                dataset = Dataset(path)
                if data['path'] is None:
                    self.datasets[variable][constituent.lower()]['path'] = path
                self.datasets[variable][constituent.lower()]['dataset'] = dataset
            except FileNotFoundError:
                raise FileNotFoundError(
                    f'no dataset found for "{variable}" ' f'"{constituent}" at "{path}"'
                )

        return dataset

    def _get_interpolation(
        self, variable: str, netcdf_variable: str, constituent: str, vertices: np.ndarray,
    ) -> np.ndarray:
        self._assert_vertices(vertices)

        xq = np.asarray([x + 360.0 if x < 0.0 else x for x in vertices[:, 0]]).flatten()
        yq = vertices[:, 1].flatten()
        dx = (self.x[-1] - self.x[0]) / len(self.x)
        xidx = np.logical_and(self.x >= np.min(xq) - 2.0 * dx, self.x <= np.max(xq) + 2.0 * dx)
        dy = (self.y[-1] - self.y[0]) / len(self.y)
        yidx = np.logical_and(self.y >= np.min(yq) - 2.0 * dy, self.y <= np.max(yq) + 2.0 * dy)
        xi, yi = np.meshgrid(self.x[xidx], self.y[yidx])
        xi = xi.flatten()
        yi = yi.flatten()
        dataset = self._get_dataset(variable, constituent)
        zi = dataset[netcdf_variable][yidx, xidx].flatten()
        values = griddata(
            (xi[~zi.mask], yi[~zi.mask]),
            zi[~zi.mask],
            (xq, yq),
            method='linear',
            fill_value=np.nan,
        )
        nan_idxs = np.where(np.isnan(values))
        values[nan_idxs] = griddata(
            (xi[~zi.mask], yi[~zi.mask]),
            zi[~zi.mask],
            (xq[nan_idxs], yq[nan_idxs]),
            method='nearest',
        )
        return values

    def _prepend_path(self, filename: str) -> str:
        if self.path is None:
            path = self.path
        elif isinstance(self.path, Path):
            path = self.path / filename
        else:
            path = f'{self.path}/{filename}'
        return path
