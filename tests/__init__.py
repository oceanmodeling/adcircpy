from difflib import Differ
from os import PathLike
from pathlib import Path
import re

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
    test_directory: PathLike, reference_directory: PathLike, skip_lines: {str: [int]} = None
):
    if not isinstance(test_directory, Path):
        test_directory = Path(test_directory)
    if not isinstance(reference_directory, Path):
        reference_directory = Path(reference_directory)
    if skip_lines is None:
        skip_lines = {}

    for reference_filename in reference_directory.iterdir():
        if reference_filename.is_dir():
            check_reference_directory(
                test_directory / reference_filename.name, reference_filename, skip_lines
            )
        else:
            test_filename = test_directory / reference_filename.name

            with open(test_filename) as test_file, open(reference_filename) as reference_file:
                test_lines = list(test_file.readlines())
                reference_lines = list(reference_file.readlines())

                diff = '\n'.join(Differ().compare(test_lines, reference_lines))
                message = f'"{test_filename}" != "{reference_filename}"\n{diff}'

                assert len(test_lines) == len(reference_lines), message

                lines_to_skip = set()
                for file_mask, line_indices in skip_lines.items():
                    if file_mask in str(test_filename) or re.match(
                        file_mask, str(test_filename)
                    ):
                        lines_to_skip.update(
                            line_index % len(test_lines) for line_index in line_indices
                        )

                for line_index in sorted(lines_to_skip, reverse=True):
                    del test_lines[line_index], reference_lines[line_index]

                assert '\n'.join(test_lines) == '\n'.join(reference_lines), message
