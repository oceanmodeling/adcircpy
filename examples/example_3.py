#! /usr/bin/env python
"""
This example recreates the Shinnecock Inlet test case with some added
improvements in order to demonstrate some of the capabilities of AdcircPy.

In contrast to example_1, this example generates input files that are separated
by a coldstart and hotstart phase.

The behaviour of this program is similar to the example_1.
"""

from datetime import datetime, timedelta
import pathlib
import shutil
import tarfile
import tempfile
import urllib.request
import warnings

import numpy as np

from adcircpy import AdcircMesh, AdcircRun, TidalForcing
from adcircpy.server.server_config import SlurmScript

PARENT = pathlib.Path(__file__).parent.absolute()
FORT14 = PARENT / "data/NetCDF_Shinnecock_Inlet/fort.14"


def main():
    # fetch shinnecock inlet test data
    if not FORT14.is_file():
        url = "https://www.dropbox.com/s/1wk91r67cacf132/"
        url += "NetCDF_shinnecock_inlet.tar.bz2?dl=1"
        g = urllib.request.urlopen(url)
        tmpfile = tempfile.NamedTemporaryFile()
        with open(tmpfile.name, 'b+w') as f:
            f.write(g.read())
        with tarfile.open(tmpfile.name, "r:bz2") as tar:
            tar.extractall(PARENT / "data/NetCDF_Shinnecock_Inlet/")

    # open mesh file
    mesh = AdcircMesh.open(FORT14, crs=4326)

    # let's generate the tau0 factor
    mesh.generate_tau0()

    # let's also add a mannings to the domain (constant for this example)
    mesh.mannings_n_at_sea_floor = np.full(mesh.values.shape, 0.025)

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
    driver = AdcircRun(
        mesh,
        start_date,
        end_date,
        spinup_time,
        tidal_forcing,
    )

    # request outputs
    driver.set_elevation_surface_output(sampling_frequency=timedelta(minutes=30))
    driver.set_velocity_surface_output(sampling_frequency=timedelta(minutes=30))

    # override defaults options
    driver.timestep = 4.0

    slurm_config = SlurmScript(
        account='example_account',
        slurm_ntasks=1000,
        run_name='ADCIRC_GAHM_GENERIC',
        partition='example_partition',
        duration=timedelta(hours=8),
        mail_type='all',
        mail_user='example@noaa.gov',
        log_filename='example_3.log',
        modules=['intel/2020 impi/2020 netcdf/4.7.2-parallel'],
        path_prefix='$HOME/adcirc/build'
    )
    driver.write(PARENT / "outputs/example_3", overwrite=True, slurm_config=slurm_config)

    # run parallel ADCIRC if binary is installed
    if shutil.which('padcirc') is not None:
        driver.run(PARENT / "outputs/example_3", overwrite=True)
    # run serial ADCIRC if binary is installed
    elif shutil.which('adcirc') is not None:
        driver.run(PARENT / "outputs/example_3", overwrite=True, nproc=1)
    # binaries are not installed, write to disk and exit
    else:
        msg = 'ADCIRC binaries were not found in PATH. ADCIRC will not run. '
        msg += "Writing files to disk..."
        warnings.warn(msg)
        driver.write(PARENT / "outputs/example_3", overwrite=True)


if __name__ == '__main__':
    main()
