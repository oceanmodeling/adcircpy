from tests import (
    check_reference_directory,
    DATA_DIRECTORY,
    OUTPUT_DIRECTORY,
    REFERENCE_DIRECTORY,
)

EXAMPLES_DIRECTORY = DATA_DIRECTORY.parent.parent / 'examples'


def test_example_1():
    exec(open(EXAMPLES_DIRECTORY / 'example_1.py').read())

    reference_directory = REFERENCE_DIRECTORY / 'test_example_1'
    output_directory = OUTPUT_DIRECTORY / 'example_1'

    check_reference_directory(
        output_directory, reference_directory, skip_lines={'fort.15': [0, -2]}
    )


def test_example_2():
    exec(open(EXAMPLES_DIRECTORY / 'example_2.py').read())

    reference_directory = REFERENCE_DIRECTORY / 'test_example_2'
    output_directory = OUTPUT_DIRECTORY / 'example_2'

    check_reference_directory(
        output_directory, reference_directory, skip_lines={'fort.15': [0, -2]}
    )


def test_example_3():
    exec(open(EXAMPLES_DIRECTORY / 'example_3.py').read())

    reference_directory = REFERENCE_DIRECTORY / 'test_example_3'
    output_directory = OUTPUT_DIRECTORY / 'example_3'

    check_reference_directory(
        output_directory, reference_directory, skip_lines={'fort.15': [0, -2]}
    )


def test_example_4():
    exec(open(EXAMPLES_DIRECTORY / 'example_4.py').read())

    reference_directory = REFERENCE_DIRECTORY / 'test_example_4'
    output_directory = OUTPUT_DIRECTORY / 'example_4'

    check_reference_directory(
        output_directory, reference_directory, skip_lines={'fort.15': [0, -2]}
    )
