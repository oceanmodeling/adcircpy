from stormevents import VortexTrack

from adcircpy.fort15 import Stations
from tests import check_reference_directory, OUTPUT_DIRECTORY, REFERENCE_DIRECTORY


def test_Stations():
    reference_directory = REFERENCE_DIRECTORY / 'test_Stations'
    output_directory = OUTPUT_DIRECTORY / 'test_Stations'

    if not output_directory.exists():
        output_directory.mkdir(parents=True, exist_ok=True)

    stations_1 = Stations.within_wind_swath(track=VortexTrack('florence2018'))
    stations_2 = Stations.within_wind_swath(
        track=VortexTrack('florence2018'), station_types=['ELEVATION']
    )
    stations_3 = Stations.within_wind_swath(
        track=VortexTrack('florence2018'),
        station_types={'ELEVATION': [8658120, 8670870, 8652587], 'VELOCITY': [8658120]},
    )

    stations_1.write(output_directory / 'stations_1.fort.15', overwrite=True)
    stations_2.write(output_directory / 'stations_2.fort.15', overwrite=True)
    stations_3.write(output_directory / 'stations_3.fort.15', overwrite=True)

    check_reference_directory(output_directory, reference_directory)
