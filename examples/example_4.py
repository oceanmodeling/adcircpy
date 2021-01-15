#! /usr/bin/env python

from datetime import datetime, timedelta
import pathlib

from adcircpy import AdcircMesh, AdcircRun, Tides
from adcircpy.forcing.waves.ww3 import WaveWatch3DataForcing
from adcircpy.forcing.winds.atmesh import AtmosphericMeshForcing
from adcircpy.server import SlurmConfig

PARENT = pathlib.Path(__file__).parent.absolute()
FORT14 = PARENT / "data/HSOFS/coldstart/fort.14"


def main():
    if not FORT14.is_file():
        raise RuntimeError('no HSOFS inputs found')

    mesh = AdcircMesh.open(FORT14, crs=4326)

    tidal_forcing = Tides()
    tidal_forcing.use_all()

    wind_forcing = AtmosphericMeshForcing(17, 3600)
    wave_forcing = WaveWatch3DataForcing(5, 3600)

    mesh.add_forcing(tidal_forcing)
    mesh.add_forcing(wind_forcing)
    mesh.add_forcing(wave_forcing)

    slurm = SlurmConfig(
        account='account',
        ntasks=1000,
        run_name='AdcircPy/examples/example_4.py',
        partition='partition',
        walltime=timedelta(hours=8),
        mail_type='all',
        mail_user='example@email.gov',
        log_filename='example_3.log',
        modules=['intel/2020', 'impi/2020', 'netcdf/4.7.2-parallel', 'esmf',
                 'metis'],
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
