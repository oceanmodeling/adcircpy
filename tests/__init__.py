from os import PathLike
from pathlib import Path

from filelock import FileLock
import pytest

from adcircpy.utilities import download_mesh

DATA_DIRECTORY = Path(__file__).parent.absolute().resolve() / 'data'
INPUT_DIRECTORY = DATA_DIRECTORY / 'input'
OUTPUT_DIRECTORY = DATA_DIRECTORY / 'output'
REFERENCE_DIRECTORY = DATA_DIRECTORY / 'reference'


@pytest.fixture
def shinnecock_mesh_directory(worker_id) -> Path:
    mesh_directory = INPUT_DIRECTORY / 'shinnecock'

    with FileLock(str(mesh_directory) + '.lock'):
        download_mesh(
            url='https://www.dropbox.com/s/1wk91r67cacf132/NetCDF_shinnecock_inlet.tar.bz2?dl=1',
            directory=mesh_directory,
        )

    return mesh_directory


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
