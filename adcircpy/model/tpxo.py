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


class TPXO:

    def __init__(self):
        """
        TODO: Change this part so that instead of using a
        _get_cached_directory() function, is either looks for the TPXO file
        pyenv_prefix = "/".join(sys.executable.split('/')[:-2])
        """
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
        vertices = self._assert_vertices(vertices)
        return self._get_interpolation(self.ha, constituent, vertices)

    def get_phase(self, constituent, vertices):
        vertices = self._assert_vertices(vertices)
        return self._get_interpolation(self.hp, constituent, vertices)

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
        vertices = np.asarray(vertices)
        msg = "vertices must be of shape M x 2"
        assert vertices.shape[1] == 2, msg
        return vertices

    @staticmethod
    def _query_yes_no(question, default="yes"):
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
                sys.stdout.write("Please respond with 'yes' or 'no' "
                                 "(or 'y' or 'n').\n")

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
        if self.__nc is None:

            url = 'ftp://ftp.oce.orst.edu/dist/tides/Global/'
            url += 'tpxo9_netcdf.tar.gz'
            msg = 'The TPXO database was not found the on system.'
            warnings.warn(msg)
            q = 'Would you like me to fetch the TPXO for you using the '
            q += 'internet? \n You may also cancel this operation and provide '
            q += 'the path to the h_tpxo9.v1.nc file using the TPXO_NCFILE '
            q += 'environment variable. '
            q += 'You may download the TPXO file from: {}'.format(url)
            a = self._query_yes_no(q)
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
        return ['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1', 'Mm', 'Mf',
                'M4', 'MN4', 'MS4', '2N2', 'S1']
