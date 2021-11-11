# ! /usr/bin/env python

from adcircpy.forcing.winds.best_track import BestTrackForcing
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
