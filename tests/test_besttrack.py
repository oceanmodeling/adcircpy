# ! /usr/bin/env python
from dateutil.parser import parse as parse_date
import pytest
import pytest_socket

from adcircpy.forcing.winds.best_track import BestTrackForcing, VortexForcing
from tests import (
    check_reference_directory,
    INPUT_DIRECTORY,
    OUTPUT_DIRECTORY,
    REFERENCE_DIRECTORY,
)


def test_from_fort22():
    input_directory = INPUT_DIRECTORY / 'test_from_fort22'
    output_directory = OUTPUT_DIRECTORY / 'test_from_fort22'
    reference_directory = REFERENCE_DIRECTORY / 'test_from_fort22'

    if not output_directory.exists():
        output_directory.mkdir(parents=True, exist_ok=True)

    best_track = BestTrackForcing.from_fort22(
        fort22=input_directory / 'irma2017_fort.22', nws=20,
    )

    assert best_track.storm_id == 'AL112017'
    assert best_track.name == 'IRMA'

    best_track.write(output_directory / 'irma2017_fort.22', overwrite=True)

    check_reference_directory(output_directory, reference_directory)


def test_from_atcf(mocker):
    input_directory = INPUT_DIRECTORY / 'test_from_atcf'
    output_directory = OUTPUT_DIRECTORY / 'test_from_atcf'
    reference_directory = REFERENCE_DIRECTORY / 'test_from_atcf'

    if not output_directory.exists():
        output_directory.mkdir(parents=True, exist_ok=True)

    best_track = BestTrackForcing.from_atcf_file(
        atcf=input_directory / 'florence2018_atcf.trk', nws=8,
    )

    assert best_track.storm_id is None
    assert best_track.name == 'WRT00001'

    best_track.write(output_directory / 'florence2018_fort.22', overwrite=True)

    check_reference_directory(output_directory, reference_directory)

    mocker.patch('matplotlib.pyplot.show')
    best_track.plot_track(show=True)


def test_recompute_velocity():
    output_directory = OUTPUT_DIRECTORY / 'test_recompute_velocity'
    reference_directory = REFERENCE_DIRECTORY / 'test_recompute_velocity'

    if not output_directory.exists():
        output_directory.mkdir(parents=True, exist_ok=True)

    best_track = BestTrackForcing('irma2017', nws=8)

    best_track.dataframe['latitude'][5] += 0.1
    best_track.dataframe['longitude'][5] -= 0.1

    best_track.write(output_directory / 'irma2017_fort.22', overwrite=True)

    check_reference_directory(output_directory, reference_directory)


def test_vortex_types():
    output_directory = OUTPUT_DIRECTORY / 'test_vortex_types'
    reference_directory = REFERENCE_DIRECTORY / 'test_vortex_types'

    if not output_directory.exists():
        output_directory.mkdir(parents=True, exist_ok=True)

    file_decks = {
        'a': {
            'start_date': parse_date('2018-09-11 06:00'),
            'end_date': None,
            'record_types': ['OFCL', 'HWRF', 'HMON', 'CARQ'],
        },
        'b': {
            'start_date': parse_date('2018-09-11 06:00'),
            'end_date': parse_date('2018-09-18 06:00'),
            'record_types': ['BEST'],
        },
    }

    for file_deck, values in file_decks.items():
        for record_type in values['record_types']:
            cyclone = VortexForcing(
                'al062018',
                start_date=values['start_date'],
                end_date=values['end_date'],
                file_deck=file_deck,
                record_type=record_type,
            )

            cyclone.write(
                output_directory / f'{file_deck}-deck_{record_type}.txt', overwrite=True
            )

    check_reference_directory(output_directory, reference_directory)


@pytest.mark.disable_socket
def test_no_internet():
    input_directory = INPUT_DIRECTORY / 'test_no_internet'
    output_directory = OUTPUT_DIRECTORY / 'test_no_internet'
    reference_directory = REFERENCE_DIRECTORY / 'test_no_internet'

    if not output_directory.exists():
        output_directory.mkdir(parents=True, exist_ok=True)

    with pytest.raises(pytest_socket.SocketBlockedError):
        VortexForcing(storm='al062018', start_date='20180911', end_date=None)

    vortex = VortexForcing.from_fort22(input_directory / 'fort.22')
    vortex.write(output_directory / 'fort.22', overwrite=True)

    check_reference_directory(output_directory, reference_directory)
