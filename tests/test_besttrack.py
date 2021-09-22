# ! /usr/bin/env python
from copy import copy

from dateutil.parser import parse as parse_date
import pytest
from pytest_socket import SocketBlockedError

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


def test_from_atcf():
    input_directory = INPUT_DIRECTORY / 'test_from_atcf'
    output_directory = OUTPUT_DIRECTORY / 'test_from_atcf'
    reference_directory = REFERENCE_DIRECTORY / 'test_from_atcf'

    if not output_directory.exists():
        output_directory.mkdir(parents=True, exist_ok=True)

    best_track = BestTrackForcing.from_atcf_file(
        atcf=input_directory / 'florence2018_atcf.trk', nws=8,
    )

    assert best_track.storm_id == 'BT02008'
    assert best_track.name == 'WRT00001'

    best_track.write(output_directory / 'florence2018_fort.22', overwrite=True)

    check_reference_directory(output_directory, reference_directory)


def test_plot_besttrack(mocker):
    input_directory = INPUT_DIRECTORY / 'test_plot_besttrack'

    best_track = BestTrackForcing.from_atcf_file(
        atcf=input_directory / 'florence2018_atcf.trk', nws=8,
    )

    mocker.patch('matplotlib.pyplot.show')
    best_track.plot_track(show=True, coastline=False)


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

    with pytest.raises((ConnectionError, SocketBlockedError)):
        VortexForcing(storm='florence2018')

    with pytest.raises((ConnectionError, SocketBlockedError)):
        VortexForcing(storm='al062018', start_date='20180911', end_date=None)

    vortex_1 = VortexForcing.from_fort22(input_directory / 'fort.22')
    vortex_1.write(output_directory / 'vortex_1.22', overwrite=True)

    vortex_2 = VortexForcing.from_fort22(vortex_1.filename)
    vortex_2.write(output_directory / 'vortex_2.22', overwrite=True)

    vortex_3 = copy(vortex_1)
    vortex_3.write(output_directory / 'vortex_3.22', overwrite=True)

    assert vortex_1 == vortex_2
    assert vortex_1 == vortex_3

    check_reference_directory(output_directory, reference_directory)
