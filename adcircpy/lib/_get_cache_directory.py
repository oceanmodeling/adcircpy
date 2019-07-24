import os


def _get_cache_directory():
    cachedir = os.getenv('LOCALAPPDATA')
    if cachedir is None:
        cachedir = os.getenv('HOME')+'/.cache/adcircpy'
    else:
        cachedir += '/adcircpy'
    os.makedirs(cachedir, exist_ok=True)
    return cachedir
