import os
import wget
import tarfile


def get_cache_dir():
    cachedir = os.getenv('LOCALAPPDATA')
    if cachedir is None:
        cachedir = os.getenv('HOME')+'/.cache/AdcircPy'
    else:
        cachedir += '/AdcircPy'
    return cachedir


def __init_TPXO_cache():
    cachedir = get_cache_dir()
    tpxo_path = cachedir + "/h_tpxo9.v1.nc"
    os.makedirs(cachedir, exist_ok=True)
    if os.path.isfile(cachedir+"/h_tpxo9.v1.nc") is False:
        print('Building TPXO database cache on {},'.format(tpxo_path)
              + ' please wait... (This will only happen the first time you run'
              + ' this software)')
        url = 'ftp://ftp.oce.orst.edu/dist/tides/Global/tpxo9_netcdf.tar.gz'
        if os.path.isfile(cachedir+"/tpxo9_netcdf.tar.gz") is False:
            wget.download(url, out=cachedir+"/tpxo9_netcdf.tar.gz")
        tpxo = tarfile.open(cachedir+"/tpxo9_netcdf.tar.gz")
        tpxo.extract('h_tpxo9.v1.nc', path=cachedir)
        tpxo.close()
        os.remove(cachedir+"/tpxo9_netcdf.tar.gz")
    return tpxo_path


def anthem():
    try:
        import pygame
    except ImportError:
        raise ImportError("Richard Stallman says: you can install pygame from "
                          + "pip if you want to hear me sing!")

    cachedir = get_cache_dir()
    if os.path.isfile(cachedir+"/free-sotfware-song.au") is False:
        wget.download('https://www.gnu.org/music/free-software-song.au',
                      out=cachedir+"/free-software-song.au",
                      bar=None)
    pygame.mixer.init()
    song = pygame.mixer.Sound(cachedir+"/free-sotfware-song.au")
    song.play()
    pygame.quit()


from AdcircPy.utils.PBS import PBS  # noqa:E402
from AdcircPy.utils.ServerConfiguration import ServerConfiguration  # noqa:E402
# from AdcircPy.utils.CTFS import CTFS  # noqa:E402

__all__ = ['PBS',
           'ServerConfiguration']

if __name__ == '__main__':
    anthem()
