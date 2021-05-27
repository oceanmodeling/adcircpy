#! /usr/bin/env python
"""
This example recreates the Shinnecock Inlet test case in it's canonical form.
The only difference between this example and the canonical case is the model
start date. The actual start date for the Shinnecock Inlet canonical example is
not documented, therefore the NCDATE that appears on the canonical files is
used.

The output files are found in outputs/example_1

If either the padcirc or adcirc binaries are found in PATH, this script will
execute the run. Otherwise it will write the files and exit.
"""

from datetime import datetime, timedelta
import pathlib
import shutil
import tarfile
import tempfile
import urllib.request
import warnings

from adcircpy import AdcircMesh, AdcircRun, Tides

PARENT = pathlib.Path(__file__).parent.absolute()
FORT14 = PARENT / 'data/NetCDF_Shinnecock_Inlet/fort.14'


def main():
    # fetch shinnecock inlet test data
    if not FORT14.is_file():
        url = 'https://www.dropbox.com/s/1wk91r67cacf132/'
        url += 'NetCDF_shinnecock_inlet.tar.bz2?dl=1'
        g = urllib.request.urlopen(url)
        tmpfile = tempfile.NamedTemporaryFile()
        with open(tmpfile.name, 'b+w') as f:
            f.write(g.read())
        with tarfile.open(tmpfile.name, 'r:bz2') as tar:
            tar.extractall(PARENT / 'data/NetCDF_Shinnecock_Inlet/')

    # open mesh file
    mesh = AdcircMesh.open(FORT14, crs=4326)

    # init tidal forcing and setup requests
    tidal_forcing = Tides()
    tidal_forcing.use_constituent('M2')
    tidal_forcing.use_constituent('N2')
    tidal_forcing.use_constituent('S2')
    tidal_forcing.use_constituent('K1')
    tidal_forcing.use_constituent('O1')

    mesh.add_forcing(tidal_forcing)

    # set simulation dates
    start_date = datetime(2015, 12, 14)
    end_date = start_date + timedelta(days=5)

    # instantiate AdcircRun object.
    driver = AdcircRun(mesh, start_date, end_date)

    # request outputs
    driver.set_elevation_surface_output(sampling_rate=timedelta(minutes=30))
    driver.set_velocity_surface_output(sampling_rate=timedelta(minutes=30))

    # override the AdcircPy defaults so that the fort.15
    # matches the original Shinnecock test case options
    driver.timestep = 6.0
    driver.DRAMP = 2.0
    driver.TOUTGE = 3.8
    driver.TOUTGV = 3.8
    driver.smagorinsky = False
    driver.horizontal_mixing_coefficient = 5.0
    driver.gwce_solution_scheme = 'semi-implicit-legacy'

    # run parallel ADCIRC if binary is installed
    if shutil.which('padcirc') is not None:
        driver.run(PARENT / 'outputs/example_1', overwrite=True)
    # run serial ADCIRC if binary is installed
    elif shutil.which('adcirc') is not None:
        driver.run(PARENT / 'outputs/example_1', overwrite=True, nproc=1)
    # binaries are not installed, write to disk and exit
    else:
        msg = 'ADCIRC binaries were not found in PATH. ADCIRC will not run. '
        msg += 'Writing files to disk...'
        warnings.warn(msg)
        driver.write(PARENT / 'outputs/example_1', overwrite=True)


if __name__ == '__main__':
    main()
