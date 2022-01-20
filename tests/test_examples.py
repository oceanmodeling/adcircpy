from tests import (
    DATA_DIRECTORY,
    OUTPUT_DIRECTORY,
    REFERENCE_DIRECTORY,
    check_reference_directory,
)

EXAMPLES_DIRECTORY = DATA_DIRECTORY.parent.parent / 'examples'


def test_example_3():
    exec(open(EXAMPLES_DIRECTORY / 'example_3.py').read())

    reference_directory = REFERENCE_DIRECTORY / 'example_3'
    output_directory = OUTPUT_DIRECTORY / 'example_3'

    check_reference_directory(
        output_directory, reference_directory, skip_lines={'fort.15': [0, -2]}
    )
