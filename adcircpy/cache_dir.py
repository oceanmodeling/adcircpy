import os


def cache_dir():
    cachedir = os.getenv('LOCALAPPDATA')
    if cachedir is None:
        cachedir = os.getenv('HOME')+'/.cache/adcircpy'
    else:
        cachedir += '/adcircpy'
    os.makedirs(cachedir, exist_ok=True)
    return cachedir
