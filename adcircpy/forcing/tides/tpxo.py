import os
import pathlib
import sys
import tarfile
import tempfile

from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import griddata
import wget


class TPXO:
    """
    Egbert, Gary D., and Svetlana Y. Erofeeva. "Efficient inverse modeling of barotropic ocean tides." Journal of Atmospheric and Oceanic Technology 19.2 (2002): 183-204  # noqa:E501
    """

    def __init__(self):

        file = os.getenv('TPXO_NCFILE')
        if file is not None:
            self._nc = Dataset(file)
            return

        else:
            prefix = os.path.dirname(os.path.dirname(sys.executable))
            file = pathlib.Path(prefix) / 'lib/h_tpxo9.v1.nc'

        if isinstance(file, pathlib.Path):
            if not file.is_file():
                self._fetch_tpxo_file(prefix, file)
            self._nc = Dataset(file)
            return

        msg = "No TPXO file found. You need to register and request a "
        msg += "copy of the TPXO9 netcdf file (specifically h_tpxo9.v1.nc)"
        msg += " from the authors at https://www.tpxo.net. Once you obtain"
        msg += " this copy, set the environment variable TPXO_NCFILE "
        msg += "to point to the path of the h_tpxo9.v1.nc file. \n"
        msg += "You may also install this file manually by placing it on "
        msg += f"the {str(file)} path."
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
                np.logical_and(x >= np.min(_x) - 2 * dx,
                               x <= np.max(_x) + 2 * dx),
                np.logical_and(y >= np.min(_y) - 2 * dy,
                               y <= np.max(_y) + 2 * dy)))
        return griddata(
            (x[_idx], y[_idx]), array[_idx], (_x, _y), method='nearest')

    def _assert_vertices(self, vertices):
        msg = "vertices must be of shape M x 2"
        assert vertices.shape[1] == 2, msg

    @staticmethod
    def _fetch_tpxo_file(prefix: str, file: str):
        url = "https://www.dropbox.com/s/uc44cbo5s2x4n93/"
        url += "h_tpxo9.v1.tar.gz?dl=1"

        def query_yes_no(question: str, default: str = "yes") -> bool:
            """
            Ask a yes/no question via raw_input() and return their answer.

            :param question: string presented to the user
            :param default: presumed answer if the user just hits <Enter>; must be "yes" (the default), "no", or None (meaning an answer is required of the user)
            :returns: whether 'yes' or 'no' was selected by the user
            """

            valid = {
                "yes": True,
                "y": True,
                "ye": True,
                "no": False,
                "n": False
            }
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
                    return valid[default]
                elif choice in valid.keys():
                    return valid[choice]
                else:
                    sys.stdout.write(
                        "Please respond with 'yes' or 'no' (or 'y' or 'n').\n")

        q = "******* PLEASE READ *******\n"
        q += "A function that is being invoked requires the TPXO file.\n"
        q += 'This software can automatically fetch the TPXO file for you usin'
        q += 'g your internet connection.\n'
        q += 'You may also cancel this operation and provide the path to the '
        q += 'h_tpxo9.v1.nc file using the TPXO_NCFILE environment variable.\n'
        q += "By downloading this file and using this software, you are "
        q += "accepting the licensing agreement for the TPXO file found here:"
        q += "\nhttps://drive.google.com/file/d/1f00WojHqu7_VE5Hg9OdiVBjymH76d"
        q += "FCy/view\n"
        q += 'If you accept the agreement, you may also download the TPXO '
        q += f"file from: {url}\n"
        q += "Would you like this software to fetch and stage the TPXO file "
        q += "from the internet now?\n"
        a = query_yes_no(q)
        if a is False:
            raise RuntimeError('No TPXO database found.')
        else:
            msg = f'Downloading TPXO database to {str(file)}'
            msg += ' please wait... \n'
            msg += 'The h_tpxo.v1.nc file will occupy about 1.2G.'
            print(msg)
            tmpdir = tempfile.TemporaryDirectory()
            _tmpdir = pathlib.Path(tmpdir.name)
            wget.download(url, out=str(_tmpdir / "h_tpxo9.v1.tar.gz"))
            with tarfile.open(_tmpdir / "h_tpxo9.v1.tar.gz") as f:
                f.extract('h_tpxo9.v1.nc', path=prefix + '/lib')

    @property
    def _nc(self):
        return self.__nc

    @_nc.setter
    def _nc(self, nc):
        self.__nc = nc


def install():
    prefix = "/".join(sys.executable.split('/')[:-2])
    file = pathlib.Path(prefix) / 'lib/h_tpxo9.v1.nc'
    TPXO._fetch_tpxo_file(prefix, file)
