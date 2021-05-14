#! /usr/bin/env python

from datetime import datetime, timedelta
import pathlib
import tarfile
import tempfile
import urllib.request

from adcircpy import AdcircMesh, AdcircRun, Tides
from adcircpy.forcing.waves.ww3 import WaveWatch3DataForcing
from adcircpy.forcing.winds.atmesh import AtmosphericMeshForcing
from adcircpy.server import SlurmConfig

PARENT = pathlib.Path(__file__).parent.absolute()
FORT14 = PARENT / "data/NetCDF_Shinnecock_Inlet/fort.14"


def main():
    if not FORT14.is_file():
        url = 'https://www.dropbox.com/s/1wk91r67cacf132/NetCDF_shinnecock_inlet.tar.bz2?dl=1'
        g = urllib.request.urlopen(url)
        tmpfile = tempfile.NamedTemporaryFile()
        with open(tmpfile.name, 'b+w') as f:
            f.write(g.read())
        with tarfile.open(tmpfile.name, "r:bz2") as tar:
            tar.extractall(PARENT / "data/NetCDF_Shinnecock_Inlet/")

    mesh = AdcircMesh.open(FORT14, crs=4326)

    tidal_forcing = Tides()
    tidal_forcing.use_all()

    wind_forcing = AtmosphericMeshForcing(
            filename='Wind_HWRF_SANDY_Nov2018_ExtendedSmoothT.nc',
            nws=17,
            interval_seconds=3600,
    )
    wave_forcing = WaveWatch3DataForcing(
            filename='ww3.HWRF.NOV2018.2012_sxy.nc',
            nrs=5,
            interval_seconds=3600,
    )

    mesh.add_forcing(tidal_forcing)
    mesh.add_forcing(wind_forcing)
    mesh.add_forcing(wave_forcing)

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
    )

    driver.write(PARENT / "outputs/example_4", overwrite=True)


if __name__ == '__main__':
    main()
