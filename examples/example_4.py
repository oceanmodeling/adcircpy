#! /usr/bin/env python

from datetime import datetime, timedelta
import pathlib
import tarfile
import tempfile
import urllib.request

from adcircpy import AdcircMesh, AdcircRun, Tides
from adcircpy.forcing.waves.base import WaveForcing
from adcircpy.forcing.winds.base import WindForcing
from adcircpy.server import SlurmConfig

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

    # init tidal forcing and setup requests
    tidal_forcing = Tides()
    tidal_forcing.use_all()

    mesh.add_forcing(tidal_forcing)

    # instantiate AdcircRun object.
    slurm = SlurmConfig(
        account='account',
        ntasks=1000,
        run_name='AdcircPy/examples/example_3.py',
        partition='partition',
        walltime=timedelta(hours=8),
        mail_type='all',
        mail_user='example@email.gov',
        log_filename='example_3.log',
        modules=['intel/2020', 'impi/2020', 'netcdf/4.7.2-parallel'],
        path_prefix='$HOME/adcirc/build'
    )
    driver = AdcircRun(
        mesh=mesh,
        start_date=datetime.now(),
        end_date=timedelta(days=7),
        spinup_time=timedelta(days=5),
        server_config=slurm,
        wind_forcing=WindForcing(17, 3600),
        wave_forcing=WaveForcing(5, 3600),
    )

    driver.write(PARENT / "outputs/example_4", overwrite=True)


if __name__ == '__main__':
    main()
