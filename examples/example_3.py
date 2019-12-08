#! /usr/bin/env python
"""
Hurrican Sandy GAHM on HSOFS
"""

import pathlib
from datetime import datetime, timedelta
from adcircpy import AdcircMesh, TidalForcing, AdcircRun

try:
    import colored_traceback
    colored_traceback.add_hook(always=True)
except ModuleNotFoundError:
    pass


def main():

    # set paths
    parent = pathlib.Path(__file__).parent
    data = parent.absolute() / "data"
    # fort14 = parent.joinpath("data/HSOFS/CubaIkeModNOMA1enoflux.grd")
    # fort13 = parent.joinpath("data/HSOFS/fort.14")

    # download HSOFS data
    # if not fort14.is_file():
    #     url = "https://www.dropbox.com/s/1wk91r67cacf132/"
    #     url += "NetCDF_shinnecock_inlet.tar.bz2?dl=1"
    #     g = urllib.request.urlopen(url)
    #     tmpfile = tempfile.NamedTemporaryFile()
    #     with open(tmpfile.name, 'b+w') as f:
    #         f.write(g.read())
    #     with tarfile.open(tmpfile.name, "r:bz2") as tar:
    #         tar.extractall(parent.joinpath("data/NetCDF_Shinnecock_Inlet/"))

    # open mesh file
    import os
    mesh = AdcircMesh.open(
        os.getenv('HSOFS'),
        crs=4326,
        fort13=os.getenv('HSOFS_FORT13')
        )

    # turn 'on' nodal attributes on both coldstart and hotstart
    mesh.set_nodal_attribute_state('mannings_n_at_sea_floor', True, True)
    mesh.set_nodal_attribute_state(
        'primitive_weighting_in_continuity_equation', True, True)
    mesh.set_nodal_attribute_state('surface_submergence_state', True, True)

    # trun 'on' nodal attributes on hotstart only
    mesh.set_nodal_attribute_hotstart_state(
        'surface_directional_effective_roughness_length', True)
    mesh.set_nodal_attribute_hotstart_state('surface_submergence_state', True)
    mesh.set_nodal_attribute_hotstart_state('surface_canopy_coefficient', True)

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

    # instantiate AdcircRun object.
    driver = AdcircRun(
        mesh,
        tidal_forcing,
        start_date=start_date,
        end_date=end_date,
        spinup_time=spinup_time
    )

    # request surface outputs
    # timedelta(0.) enable maximum surface while disabling surface timeseries
    driver.set_elevation_surface_output(timedelta(0.))

    # override defaults to match original test case options
    driver.timestep = 2.
    driver.spinup_factor = 2./3.
    driver.smagorinsky = False
    driver.horizontal_mixing_coefficient = 10.

    # write files to disk
    outdir = data.absolute() / 'example_3'
    # driver.adcirc(outdir, overwrite=True)
    driver.padcirc(outdir, overwrite=True)


if __name__ == '__main__':
    main()
