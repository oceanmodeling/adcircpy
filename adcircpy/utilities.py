import logging
import os
from os import PathLike
from pathlib import Path
import tarfile

import wget


def download_mesh(url: str, directory: PathLike, overwrite: bool = False):
    if not isinstance(directory, Path):
        directory = Path(directory)
    if not directory.exists():
        directory.mkdir(parents=True, exist_ok=True)

    if not (directory / 'fort.14').exists() or overwrite:
        logging.info(f'downloading mesh files to {directory}')
        extract_download(url, directory, ['fort.13', 'fort.14'])

    return directory


def extract_download(
    url: str, directory: PathLike, filenames: [str] = None, overwrite: bool = False
):
    if not isinstance(directory, Path):
        directory = Path(directory)

    if filenames is None:
        filenames = []

    if not directory.exists():
        directory.mkdir(parents=True, exist_ok=True)

    temporary_filename = directory / 'temp.tar.gz'
    logging.debug(f'downloading {url} -> {temporary_filename}')
    wget.download(url, f'{temporary_filename}')
    logging.debug(f'extracting {temporary_filename} -> {directory}')
    with tarfile.open(temporary_filename) as local_file:
        if len(filenames) > 0:
            for filename in filenames:
                if filename in local_file.getnames():
                    path = directory / filename
                    if not path.exists() or overwrite:
                        if path.exists():
                            os.remove(path)
                        local_file.extract(filename, directory)
        else:
            local_file.extractall(directory)

    os.remove(temporary_filename)
