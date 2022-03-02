import logging
import os
from os import PathLike
from pathlib import Path
import sys
import tarfile

import pooch


def download_mesh(
    url: str, directory: PathLike, known_hash: str = None, overwrite: bool = False
):
    if not isinstance(directory, Path):
        directory = Path(directory)
    if not directory.exists():
        directory.mkdir(parents=True, exist_ok=True)

    if not (directory / 'fort.14').exists() or overwrite:
        logging.info(f'downloading mesh files to {directory}')
        extract_download(
            url, directory, ['fort.13', 'fort.14'], known_hash=known_hash, overwrite=overwrite
        )

    return directory


def extract_download(
    url: str,
    directory: PathLike,
    filenames: [str] = None,
    known_hash: str = None,
    overwrite: bool = False,
):
    if not isinstance(directory, Path):
        directory = Path(directory)

    if filenames is None:
        filenames = []

    if not directory.exists():
        directory.mkdir(parents=True, exist_ok=True)

    temporary_filename = directory / 'temp.tar.gz'
    logging.debug(f'downloading {url} -> {temporary_filename}')
    temporary_filename = pooch.retrieve(url, known_hash=known_hash, fname=temporary_filename)
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


def get_logger(
    name: str,
    log_filename: PathLike = None,
    file_level: int = None,
    console_level: int = None,
    log_format: str = None,
) -> logging.Logger:
    if file_level is None:
        file_level = logging.DEBUG
    if console_level is None:
        console_level = logging.INFO
    logger = logging.getLogger(name)

    # check if logger is already configured
    if logger.level == logging.NOTSET and len(logger.handlers) == 0:
        # check if logger has a parent
        if '.' in name:
            if isinstance(logger.parent, logging.RootLogger):
                for existing_console_handler in [
                    handler
                    for handler in logger.parent.handlers
                    if not isinstance(handler, logging.FileHandler)
                ]:
                    logger.parent.removeHandler(existing_console_handler)
            logger.parent = get_logger(name.rsplit('.', 1)[0])
        else:
            # otherwise create a new split-console logger
            if console_level != logging.NOTSET:
                for existing_console_handler in [
                    handler
                    for handler in logger.handlers
                    if not isinstance(handler, logging.FileHandler)
                ]:
                    logger.removeHandler(existing_console_handler)

                console_output = logging.StreamHandler(sys.stdout)
                console_output.setLevel(console_level)
                logger.addHandler(console_output)

    if log_filename is not None:
        file_handler = logging.FileHandler(log_filename)
        file_handler.setLevel(file_level)
        for existing_file_handler in [
            handler for handler in logger.handlers if isinstance(handler, logging.FileHandler)
        ]:
            logger.removeHandler(existing_file_handler)
        logger.addHandler(file_handler)

    if log_format is None:
        log_format = '[%(asctime)s] %(name)-15s %(levelname)-8s: %(message)s'
    log_formatter = logging.Formatter(log_format)
    for handler in logger.handlers:
        handler.setFormatter(log_formatter)

    return logger
