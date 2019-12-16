#! /usr/bin/env python
"""
Hurrican Sandy GAHM on HSOFS
"""
import os
import pathlib
import warnings
from datetime import datetime, timedelta
from adcircpy import AdcircMesh, TidalForcing, BestTrackForcing, AdcircRun

try:
    import colored_traceback
    colored_traceback.add_hook(always=True)
except ModuleNotFoundError:
    pass


def main():

    fort14 = os.getenv('HSOFS')
    fort13 = os.getenv('HSOFS_FORT13')
    stations = os.getenv('STATIONS')

    if fort14 is None:
        msg = 'Please set the enviroment variable\n'
        msg += 'export HSOFS=/path/to/HSOFS'
        raise IOError(msg)

    if fort13 is None:
        msg = 'Please set the enviroment variable\n'
        msg += 'export HSOFS_FORT13=/path/to/HSOFS_FORT13'
        raise IOError(msg)

    if stations is None:
        msg = 'No STATIONS environment variable was found. '
        msg += 'ADCIRC input files will be created but they will not contain '
        msg += 'output station data. Please set the STATIONS envrionment '
        msg += 'variable with station data in the same format as a fort.15. '
        msg += 'You may simple pass an existing fort.15 as it is and the '
        msg += 'station data will be extracted from it.'
        warnings.warn(msg)

    # open mesh file
    mesh = AdcircMesh.open(
        fort14,
        crs=4326,
        fort13=fort13
        )

    # turn 'on' nodal attributes on both coldstart and hotstart
    mesh.set_nodal_attribute_state('mannings_n_at_sea_floor', True, True)
    mesh.set_nodal_attribute_state(
        'primitive_weighting_in_continuity_equation', True, True)
    mesh.set_nodal_attribute_state('surface_submergence_state', True, True)

    # trun 'on' nodal attributes on hotstart only
    mesh.set_nodal_attribute_hotstart_state('surface_canopy_coefficient', True)
    mesh.set_nodal_attribute_hotstart_state(
        'surface_directional_effective_roughness_length', True)

    # init tidal forcing and add constituents
    tidal_forcing = TidalForcing()
    tidal_forcing.use_constituent('K1')
    tidal_forcing.use_constituent('O1')
    tidal_forcing.use_constituent('P1')
    tidal_forcing.use_constituent('Q1')
    tidal_forcing.use_constituent('N2')
    tidal_forcing.use_constituent('M2')
    tidal_forcing.use_constituent('S2')
    tidal_forcing.use_constituent('K2')
    tidal_forcing.use_constituent('Mf')
    tidal_forcing.use_constituent('Mm')
    tidal_forcing.use_constituent('M4')
    tidal_forcing.use_constituent('MS4')
    tidal_forcing.use_constituent('MN4')

    # set simulation dates
    spinup_time = timedelta(days=15.)
    start_date = datetime(2012, 10, 11) + spinup_time
    end_date = (start_date - spinup_time) + timedelta(days=19.25)

    # set winds
    wind_forcing = BestTrackForcing('AL182012')  # Hurricane Sandy

    # instantiate AdcircRun object.
    driver = AdcircRun(
        mesh,
        tidal_forcing,
        wind_forcing,
        start_date=start_date,
        end_date=end_date,
        spinup_time=spinup_time
    )

    # Request surface outputs. timedelta(0) enables maximum surface output
    # (e.g. maxele.63.nc) while disabling surface timeseries.
    driver.set_elevation_surface_output(timedelta(0.))

    if stations is not None:
        driver.import_stations(stations)
        driver.set_elevation_stations_output(timedelta(minutes=6))

    # override defaults to match original test case options
    driver.timestep = 2.
    driver.spinup_factor = 2./3.
    driver.smagorinsky = False
    driver.horizontal_mixing_coefficient = 10.

    # write files to disk
    parent = pathlib.Path(__file__).parent.absolute()
    driver.dump(parent / 'data/example_3', overwrite=True)


if __name__ == '__main__':
    main()
