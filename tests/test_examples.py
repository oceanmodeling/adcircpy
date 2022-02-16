import shutil
import sys

import pytest

from tests import (
    check_reference_directory,
    DATA_DIRECTORY,
    OUTPUT_DIRECTORY,
    REFERENCE_DIRECTORY,
)

EXAMPLES_DIRECTORY = DATA_DIRECTORY.parent.parent / 'examples'


@pytest.mark.skip
def test_example_1():
    reference_directory = REFERENCE_DIRECTORY / 'example_1'
    output_directory = OUTPUT_DIRECTORY / 'example_1'

    if output_directory.exists():
        shutil.rmtree(output_directory)
    output_directory.mkdir(exist_ok=True, parents=True)

    exec(open(EXAMPLES_DIRECTORY / 'example_1.py').read())

    check_reference_directory(
        output_directory, reference_directory, skip_lines={'fort.15': [0, -2]}
    )


@pytest.mark.skip
def test_example_2():
    reference_directory = REFERENCE_DIRECTORY / 'example_2'
    output_directory = OUTPUT_DIRECTORY / 'example_2'

    if output_directory.exists():
        shutil.rmtree(output_directory)
    output_directory.mkdir(exist_ok=True, parents=True)

    exec(open(EXAMPLES_DIRECTORY / 'example_2.py').read())

    check_reference_directory(
        output_directory, reference_directory, skip_lines={'fort.15': [0, -2]}
    )


@pytest.mark.skipif(sys.version_info < (3, 7), reason='test fails on Python < 3.7')
def test_example_3():
    reference_directory = REFERENCE_DIRECTORY / 'example_3'
    output_directory = OUTPUT_DIRECTORY / 'example_3'

    if output_directory.exists():
        shutil.rmtree(output_directory)
    output_directory.mkdir(exist_ok=True, parents=True)

    exec(open(EXAMPLES_DIRECTORY / 'example_3.py').read())

    check_reference_directory(
        output_directory, reference_directory, skip_lines={'fort.15': [0, -2]}
    )


@pytest.mark.skip
def test_example_4():
    reference_directory = REFERENCE_DIRECTORY / 'example_4'
    output_directory = OUTPUT_DIRECTORY / 'example_4'

    if output_directory.exists():
        shutil.rmtree(output_directory)
    output_directory.mkdir(exist_ok=True, parents=True)

    exec(open(EXAMPLES_DIRECTORY / 'example_4.py').read())

    check_reference_directory(
        output_directory,
        reference_directory,
        skip_lines={'fort.15.coldstart': [0, -2], 'fort.15.hotstart': [0, -3]},
    )
