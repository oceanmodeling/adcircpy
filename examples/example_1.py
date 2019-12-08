#! /usr/bin/env python
"""
Shinnecock Inlet example analogous to canonical test case
"""

import pathlib
import tarfile
import tempfile
import urllib.request
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
    fort14 = parent.joinpath("data/NetCDF_Shinnecock_Inlet/fort.14")

    # download_shinnecock_data
    if not fort14.is_file():
        url = "https://www.dropbox.com/s/1wk91r67cacf132/"
        url += "NetCDF_shinnecock_inlet.tar.bz2?dl=1"
        g = urllib.request.urlopen(url)
        tmpfile = tempfile.NamedTemporaryFile()
        with open(tmpfile.name, 'b+w') as f:
            f.write(g.read())
        with tarfile.open(tmpfile.name, "r:bz2") as tar:
            tar.extractall(parent.joinpath("data/NetCDF_Shinnecock_Inlet/"))

    # open mesh file
    mesh = AdcircMesh.open(fort14, crs=4326)

    # init tidal forcing and setup requests
    tidal_forcing = TidalForcing()
    tidal_forcing.use_constituent('M2')
    tidal_forcing.use_constituent('N2')
    tidal_forcing.use_constituent('S2')
    tidal_forcing.use_constituent('K1')
    tidal_forcing.use_constituent('O1')

    # set simulation dates
    start_date = datetime(2015, 12, 14)
    end_date = start_date + timedelta(days=5)

    # instantiate AdcircRun object.
    driver = AdcircRun(
        mesh,
        tidal_forcing,
        start_date=start_date,
        end_date=end_date,
    )

    # request outputs
    driver.set_elevation_surface_output(
        sampling_frequency=timedelta(minutes=30),
        spinup=True)
    driver.set_velocity_surface_output(
        sampling_frequency=timedelta(minutes=30),
        spinup=True)

    # override defaults to match original test case options
    driver.timestep = 6.
    driver.DRAMP = 2.
    driver.smagorinsky = False
    driver.horizontal_mixing_coefficient = 5.
    driver.TOUTGE = 3.8
    driver.TOUTGV = 3.8

    # write files to disk
    outdir = parent / pathlib.Path('example_1')
    driver.adcirc(outdir, overwrite=True)


if __name__ == '__main__':
    main()
