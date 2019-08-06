from scipy.interpolate import RectBivariateSpline, griddata
import numpy as np
import warnings
import os
import sys
import wget
import tarfile
from netCDF4 import Dataset
from pathlib import Path
from adcircpy.lib._get_cache_directory import _get_cache_directory


class _TPXO(object):

    __tpxo_constituents = ['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1',
                           'Mm', 'Mf', 'M4', 'MN4', 'MS4', '2N2', 'S1']

    def __init__(self):
        try:
            self.__nc = Dataset(os.getenv('TPXO_NCFILE'))
        except FileNotFoundError:
            try:
                path = str(
                    Path(_get_cache_directory() + '/' + 'h_tpxo9.v1.nc'))
                self.__nc = Dataset(path)
            except FileNotFoundError:
                self.__nc = None

    def __call__(self, constituent, vertices):
        amp = self.get_amplitude(constituent, vertices)
        phase = self.get_phase(constituent, vertices)
        return (amp, phase)

    def get_amplitude(self, constituent, vertices, method='nearest'):
        return self.__get_interp('ha', constituent, vertices, method)

    def get_phase(self, constituent, vertices, method='nearest'):
        return self.__get_interp('hp', constituent, vertices, method)

    def __get_interp(
        self,
        interp_type,
        constituent,
        vertices,
        method='nearest'
    ):
        vertices = np.asarray(vertices)
        msg = "vertices must be of shape M x 2"
        assert vertices.shape[1] == 2, msg
        if method == 'spline':
            return self.__get_spline_interp(interp_type, constituent, vertices)
        else:
            return self.__get_griddata_interp(
                interp_type, constituent, vertices, method=method)

    def __get_spline_interp(self, interp_type, constituent, vertices):
        raise NotImplementedError
        # return RectBivariateSpline(x, y, val).ev(_x, _y)

    def __get_griddata_interp(
        self,
        interp_type,
        constituent,
        vertices,
        method=None
    ):
        x = self.nc['lon_z'][:][:].reshape(self.nc['lon_z'].size)
        y = self.nc['lat_z'][:][:].reshape(self.nc['lat_z'].size)
        idx = self.tpxo_constituents.index(constituent)
        h = self.nc[interp_type][idx, :, :].reshape(
            self.nc[interp_type][idx, :, :].size)
        _x = vertices[:, 0]
        _y = vertices[:, 1]
        _x = np.asarray([_ + 360. for _ in _x if _ < 0])
        _idx = np.where(np.logical_and(
                            np.logical_and(np.min(_x) <= x, np.max(_x) >= x),
                            np.logical_and(np.min(_y) <= y, np.max(_y) >= y)))
        h = griddata((x[_idx], y[_idx]), h[_idx], (_x, _y), method=method)
        return h.data

    @property
    def nc(self):
        if self.__nc is None:
            def query_yes_no(question, default="yes"):
                """Ask a yes/no question via raw_input() and return their answer.
                "question" is a string that is presented to the user.
                "default" is the presumed answer if the user just hits <Enter>.
                    It must be "yes" (the default), "no" or None (meaning
                    an answer is required of the user).

                The "answer" return value is one of "yes" or "no".
                """
                valid = {"yes": True,   "y": True,  "ye": True,
                         "no": False,     "n": False}
                if default is None:
                    prompt = " [y/n] "
                elif default == "yes":
                    prompt = " [Y/n] "
                elif default == "no":
                    prompt = " [y/N] "
                else:
                    raise ValueError("invalid default answer: '%s'" % default)

                while 1:
                    sys.stdout.write(question + prompt)
                    choice = input().lower()
                    if default is not None and choice == '':
                        return default
                    elif choice in valid.keys():
                        return valid[choice]
                    else:
                        sys.stdout.write("Please respond with 'yes' or 'no' " \
                                         "(or 'y' or 'n').\n")
            url = 'ftp://ftp.oce.orst.edu/dist/tides/Global/'
            url += 'tpxo9_netcdf.tar.gz'
            msg = 'The TPXO database was not found the on system.'
            warnings.warn(msg)
            q = 'Would you like me to fetch the TPXO for you using the '
            q += 'internet? \n You may also cancel this operation and provide '
            q += 'the path to the h_tpxo9.v1.nc file using the TPXO_NCFILE '
            q += 'environment variable. '
            q += 'You may download the TPXO file from: {}'.format(url)
            a = query_yes_no(q)
            if a is False:
                raise RuntimeError('No TPXO database found.')
            else:
                cachedir = str(Path(_get_cache_directory()))
                tpxo_path = str(Path(cachedir + '/' + 'h_tpxo9.v1.nc'))
                os.makedirs(cachedir, exist_ok=True)
                if os.path.isfile(cachedir+"/h_tpxo9.v1.nc") is False:
                    msg = 'Downloding TPXO database to {},'.format(tpxo_path)
                    msg += ' please wait... \n'
                    msg += 'The h_tpxo.v1.nc file will occupy about 1.1G.'
                    print(msg)
                    if os.path.isfile(str(Path(
                            cachedir+"/tpxo9_netcdf.tar.gz"))) is False:
                        wget.download(url, out=cachedir+"/tpxo9_netcdf.tar.gz")
                    tpxo = tarfile.open(cachedir+"/tpxo9_netcdf.tar.gz")
                    tpxo.extract('h_tpxo9.v1.nc', path=cachedir)
                    tpxo.close()
                    os.remove(cachedir+"/tpxo9_netcdf.tar.gz")
            self.__nc = Dataset(tpxo_path)
        return self.__nc

    @property
    def tpxo_constituents(self):
        return self.__tpxo_constituents


#     def _set_TPXO_interp(self, AdcircMesh, method='nearest'):
#         """
#         Method to generate the TPXO interpolation on the open boundaries.
#         """
#         assert isinstance(AdcircMesh, Model.AdcircMesh)
#         if len(AdcircMesh.OceanBoundaries) == 0:
#             raise RuntimeError('Cannot generate TPXO parameters for mesh '
#                                'without defined ocean boundaries.')
#         x = self.nc['lon_z'][:].reshape(self.nc['lon_z'].size)
#         y = self.nc['lat_z'][:].reshape(self.nc['lat_z'].size)
#         data = list()
#         for ocean_boundary in AdcircMesh.OceanBoundaries:
#             SpatialReference = ocean_boundary['SpatialReference']
#             if not SpatialReference.IsGeographic():
#                 raise NotImplementedError(
#                     'Tidal Run only supported for meshes in '
#                     + ' geographic coordinates.')
#             _x = ocean_boundary['vertices'][:, 0]
#             _x = np.asarray([_ + 360. for _ in _x if _ < 0])
#             _y = ocean_boundary['vertices'][:, 1]
#             _idx = np.where(np.logical_and(
#                             np.logical_and(np.min(_x) <= x, np.max(_x) >= x),
#                             np.logical_and(np.min(_y) <= y, np.max(_y) >= y)))
#             constituents = CaseInsensitiveDict()
#             for constituent in self.constituents:
#                 if constituent in self.tpxo_constituents.keys():
#                     constituents[constituent] = dict()
#                     idx = self.tpxo_constituents[constituent]
#                     ha = self.nc['ha'][idx, :, :].reshape(
#                         self.nc['ha'][idx, :, :].size)
#                     hp = self.nc['hp'][idx, :, :].reshape(
#                         self.nc['hp'][idx, :, :].size)
#                     ha = griddata((x[_idx], y[_idx]), ha[_idx], (_x, _y),
#                                   method=method)
#                     hp = griddata((x[_idx], y[_idx]), hp[_idx], (_x, _y),
#                                   method=method)
#                     constituents[constituent] = np.vstack([ha, hp]).T
#             data.append(constituents)
# self.__TPXO_interp = data