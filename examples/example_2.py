#! /usr/bin/env python
"""
Shinnecock Inlet example with coldstart and hotstart parts
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
    spinup_time = timedelta(days=2)
    start_date = datetime(2015, 12, 14) + spinup_time
    end_date = start_date + timedelta(days=3)

    # instantiate AdcircRun object.
    adcirc_run = AdcircRun(
        mesh,
        tidal_forcing,
        start_date=start_date,
        end_date=end_date,
        spinup_time=spinup_time
    )

    # request outputs
    adcirc_run.set_elevation_surface_output(
        sampling_frequency=timedelta(minutes=30))
    adcirc_run.set_velocity_surface_output(
        sampling_frequency=timedelta(minutes=30))

    # override defaults options
    adcirc_run.timestep = 6.0

    outdir = parent / pathlib.Path('example_2')
    outdir.mkdir(parents=True, exist_ok=True)
    adcirc_run.dump(outdir, overwrite=True)


if __name__ == '__main__':
    main()
