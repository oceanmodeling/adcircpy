import os
import pathlib
import sys

from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import griddata

# https://icdc.cen.uni-hamburg.de/en/hamtide.html
base_url = 'https://icdc.cen.uni-hamburg.de/thredds/dodsC/ftpthredds/hamtide/'


class HamtideResource:

    def __set__(self, obj, resource):
        if resource is None:
            resource = {'elevation': {}, 'velocity': {}}
            for constituent in obj.constituents:
                for key in resource.keys():
                    resource[key].update({constituent: None})
        else:
            raise NotImplementedError('Check that static files exist.')
            _resource = pathlib.Path(resource)
            for file in _resource.glob('*.nc'):
                print(file)

        obj.__dict__['resource'] = resource

    def __get__(self, obj, val):
        return obj.__dict__['resource']


class HAMTIDE:
    '''
    Taguchi, E., D. Stammer and W. Zahel (2010), Estimation of deep ocean tidal energy dissipation based on the high-resolution data-assimilative HAMTIDE model (to be submitted to J. Geophys. Res.). 
    '''

    _resource = HamtideResource()

    def __init__(self, resource=None):
        self._resource = resource

    def __call__(self, constituent, vertices):
        amp = self.get_amplitude(constituent, vertices)
        phase = self.get_phase(constituent, vertices)
        return (amp, phase)

    def get_amplitude(self, constituent, vertices):
        vertices = np.asarray(vertices)
        self._assert_vertices(vertices)
        return self._get_interpolation(
            'elevation', 'AMPL', constituent, vertices)

    def get_phase(self, constituent, vertices):
        vertices = np.asarray(vertices)
        self._assert_vertices(vertices)
        return self._get_interpolation(
            'elevation', 'PHAS', constituent, vertices)

    @property
    def x(self):
        if not hasattr(self, '_x'):
            self._x = Dataset(base_url + 'k2.hamtide11a.nc')['LON'][:].data
        return self._x

    @property
    def y(self):
        if not hasattr(self, '_y'):
            self._y = Dataset(base_url + 'k2.hamtide11a.nc')['LAT'][:].data
        return self._y

    @property
    def constituents(self):
        return ['S2', 'Q1', 'P1', 'O1', 'N2', 'M2', 'K2', 'K1']

    def _get_resource(self, variable, constituent) -> Dataset:
        resource = self._resource[variable][constituent]
        if resource is not None:
            return Dataset(resource)
        if variable == 'elevation':
            fname = f'{constituent.lower()}.hamtide11a.nc'
        if variable == 'velocity':
            fname = f'HAMcurrent11a_{constituent.lower()}.nc'
        return Dataset(base_url + fname)

    def _get_interpolation(self, phys_var, ncvar, constituent, vertices):
        xq = np.asarray(
            [x + 360. if x < 0. else x for x in vertices[:, 0]]).flatten()
        yq = vertices[:, 1].flatten()

        xidx = np.logical_and(self.x >= np.min(xq), self.x <= np.max(xq))
        yidx = np.logical_and(self.y >= np.min(yq), self.y <= np.max(yq))

        xi, yi = np.meshgrid(self.x[xidx], self.y[yidx])
        xi = xi.flatten()
        yi = yi.flatten()
        zi = self._get_resource(
            phys_var, constituent)[ncvar][yidx, xidx].flatten()

        return griddata((xi[~zi.mask], yi[~zi.mask]), zi[~zi.mask], (xq, yq),
                        method='nearest') * 0.01

    def _assert_vertices(self, vertices):
        msg = "vertices must be of shape M x 2"
        assert vertices.shape[1] == 2, msg
