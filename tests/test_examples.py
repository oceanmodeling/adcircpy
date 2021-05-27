from pathlib import Path
import unittest

from examples.example_3 import main, PARENT


class TestExamples(unittest.TestCase):
    def test_example_3(self):
        main()

        reference_fort15_coldstart = (
            Path(__file__) / '../reference/test_examples/example_3/fort.15.coldstart'
        )
        output_fort15_coldstart = PARENT / 'outputs/example_3/fort.15.coldstart'

        reference_fort15_hotstart = (
            Path(__file__) / '../reference/test_examples/example_3/fort.15.hotstart'
        )
        output_fort15_hotstart = PARENT / 'outputs/example_3/fort.15.hotstart'

        with open(output_fort15_coldstart) as output_fort15, open(
            reference_fort15_coldstart
        ) as reference_fort15:
            output_lines = output_fort15.readlines()
            reference_lines = reference_fort15.readlines()
            self.assertEqual(reference_lines[1:-1], output_lines[1:-1])

        with open(output_fort15_hotstart) as output_fort15, open(
            reference_fort15_hotstart
        ) as reference_fort15:
            output_lines = output_fort15.readlines()
            reference_lines = reference_fort15.readlines()
            self.assertEqual(reference_lines[1:-1], output_lines[1:-1])


if __name__ == '__main__':
    unittest.main()
