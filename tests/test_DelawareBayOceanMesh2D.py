#!/usr/bin/env python
import os
from datetime import timedelta
from pathlib import Path
import unittest
from adcircpy.mesh import AdcircMesh
from adcircpy.model import BestTrackForcing, TidalForcing, AdcircRun

DelawareBayMesh = str(Path(os.getenv('DelawareBayOceanMesh2D_fort14')))
HsofsSandyFort15 = str(Path(os.getenv('HSOFS_SANDY_HOTSTART')))


class DelawareBayOceanMesh2DTestCase(unittest.TestCase):

    def test_Isabel2003BestTrack(self):

        #  ----------- instantiante run
        adcirc_run = AdcircRun()
        adcirc_run.AdcircMesh = AdcircMesh.open(DelawareBayMesh, 4326)
        adcirc_run.AdcircMesh.import_nodal_attributes(
            os.getenv('DelawareBayOceanMesh2D_fort13'))
        for attribute in adcirc_run.AdcircMesh.get_nodal_attribute_names():
            adcirc_run.AdcircMesh.set_nodal_attribute_state(
                attribute, True, True)

        #  ----------- add tidal forcing
        adcirc_run.TidalForcing = TidalForcing()
        adcirc_run.TidalForcing.use_major()

        #  ----------- add wind forcing
        adcirc_run.WindForcing = BestTrackForcing()
        adcirc_run.WindForcing.storm_id = 'AL132003'
        adcirc_run.WindForcing.remove_TS()
        # adcirc_run.WindForcing.remove_EX()

        #  ----------- request outputs
        adcirc_run.copy_fort15_stations(HsofsSandyFort15)
        adcirc_run.set_elevation_stations_output(timedelta(minutes=6.))
        adcirc_run.set_elevation_global_output(timedelta(minutes=15.))

        #  ----------- set additional run options
        adcirc_run.start_date = adcirc_run.WindForcing.start_date
        adcirc_run.end_date = adcirc_run.WindForcing.end_date
        adcirc_run.spinup_time = timedelta(days=7.)
        adcirc_run.DTDP = 2.
        adcirc_run.ESLM = -0.05

        #  ----------- dump run
        output_directory = str(Path('test_DelawareBayOceanMesh2D'))
        os.makedirs(output_directory, exist_ok=True)
        adcirc_run.dump(output_directory, overwrite=True)


if __name__ == '__main__':
    unittest.main()

