# ! /usr/bin/env python

from adcircpy.forcing.winds.best_track import BestTrackForcing
from tests import INPUT_DIRECTORY


def test_plot_besttrack(mocker):
    input_directory = INPUT_DIRECTORY / 'test_plot_besttrack'

    best_track = BestTrackForcing.from_atcf_file(
        atcf=input_directory / 'florence2018_atcf.trk', nws=8,
    )

    mocker.patch('matplotlib.pyplot.show')
    best_track.plot_track(coastline=False)
