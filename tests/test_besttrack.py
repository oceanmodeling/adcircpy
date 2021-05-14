from pathlib import Path
import unittest

from adcircpy.forcing.winds import BestTrackForcing

DATA_DIRECTORY = Path(__file__).parent / 'data'
INPUT_DIRECTORY = DATA_DIRECTORY / 'input'
OUTPUT_DIRECTORY = DATA_DIRECTORY / 'output'
REFERENCE_DIRECTORY = DATA_DIRECTORY / 'reference'


class TestBestTrack(unittest.TestCase):
    def test_from_fort22(self):
        input_filename = INPUT_DIRECTORY / 'test_besttrack' / 'irma2017_fort.22'
        output_filename = OUTPUT_DIRECTORY / 'test_besttrack' / 'irma2017_fort.22'
        reference_filename = REFERENCE_DIRECTORY / 'test_besttrack' / 'irma2017_fort.22'

        if not output_filename.parent.exists():
            output_filename.parent.mkdir(parents=True, exist_ok=True)

        best_track = BestTrackForcing.from_fort22(
                fort22=input_filename,
                nws=20,
        )

        assert best_track.storm_id == 'AL112017'
        assert best_track.name == 'IRMA'

        best_track.write(output_filename, overwrite=True)

        with open(output_filename) as test_file:
            with open(reference_filename) as reference_file:
                assert test_file.read() == reference_file.read()
