import logging
import os
from os import PathLike
from pathlib import Path
import tarfile

from filelock import FileLock
import pytest
import wget

DATA_DIRECTORY = Path(__file__).parent.absolute().resolve() / 'data'
INPUT_DIRECTORY = DATA_DIRECTORY / 'input'
OUTPUT_DIRECTORY = DATA_DIRECTORY / 'output'
REFERENCE_DIRECTORY = DATA_DIRECTORY / 'reference'


@pytest.fixture
def shinnecock_mesh_directory(worker_id) -> Path:
    mesh_directory = INPUT_DIRECTORY / 'shinnecock'

    download_mesh(
        url='https://www.dropbox.com/s/1wk91r67cacf132/NetCDF_shinnecock_inlet.tar.bz2?dl=1',
        directory=mesh_directory,
    )

    return mesh_directory


def download_mesh(url: str, directory: PathLike, overwrite: bool = False):
    if not isinstance(directory, Path):
        directory = Path(directory)
    if not directory.exists():
        directory.mkdir(parents=True, exist_ok=True)

    if not (directory / 'fort.14').exists() or overwrite:
        with FileLock(str(directory) + '.lock'):
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


def check_reference_directory(
    test_directory: PathLike, reference_directory: PathLike, skip_lines: int = None
):
    if not isinstance(test_directory, Path):
        test_directory = Path(test_directory)
    if not isinstance(reference_directory, Path):
        reference_directory = Path(reference_directory)
    if skip_lines is None:
        skip_lines = 0

    for reference_filename in reference_directory.iterdir():
        if reference_filename.is_dir():
            check_reference_directory(
                test_directory / reference_filename.name, reference_filename, skip_lines
            )
        else:
            test_filename = test_directory / reference_filename.name
            with open(test_filename) as test_file, open(reference_filename) as reference_file:
                assert (
                    test_file.readlines()[skip_lines:]
                    == reference_file.readlines()[skip_lines:]
                ), f'"{test_filename}" != "{reference_filename}"'
