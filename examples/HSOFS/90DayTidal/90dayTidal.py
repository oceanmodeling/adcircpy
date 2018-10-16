#! /usr/bin/env python
import os
from datetime import datetime, timedelta
from AdcircPy import AdcircPy
from AdcircPy import ElevationGlobalOutput as EGO
from AdcircPy import ElevationStationsOutput as ESO

MESH_PATH='./fort.14'
FORT13_PATH='./fort.13'
FORT15_PATH='./fort.15'

class _90DayTidal(object):
  """
  Example class for generating the 90 day Tidal run.
  The outputs provided in the example folder correspond to the HSOFS mesh.
  The class looks for the mesh path specified on the environment variable
  FORT14_PATH. Optionally, it imports the elevation stations from a fort.15
  containing the NOUTE variable specifying elevation stations.
  If a fort.15 is not provided, it will still create the forcing file,
  but will not include station outputs.
  """
  
  def __init__(self):
    self._init_mesh()
    self._init_dates()
    self._init_elevation_station_output_request()
    self._init_elevation_global_output_request()
    self._init_TidalRun()
    self._dump_TidalRun()

  def _init_mesh(self):
    self.Mesh = AdcircPy.read_mesh(fort14=MESH_PATH,
                                    fort13=FORT13_PATH)

  def _init_dates(self):
    self.spinup_date = datetime(2013, 8, 1, 0, 0, 0)
    self.start_time  = self.spinup_date + timedelta(days=30)
    self.end_time    = self.spinup_date + timedelta(days=90)
  
  def _init_elevation_station_output_request(self):
    self.ElevationStationsOutput = ESO.from_fort15(FORT15_PATH,
                                       harmonic_analysis=True,
                                       spinup=True)

  def _init_elevation_global_output_request(self):
    self.ElevationGlobalOutput = EGO(harmonic_analysis=True,
                                     spinup=True)

  def _init_TidalRun(self):
    self.TidalRun = self.Mesh.TidalRun(self.start_time, self.end_time,
                                        spinup_date=self.spinup_date,
                                        ESO=self.ElevationStationsOutput,
                                        EGO=self.ElevationGlobalOutput,
                                        constituents=['K1', 'O1', 'P1', 'Q1',
                                                      'N2', 'M2', 'S2', 'K2',
                                                      'Mf', 'Mm', 'M4', 'MS4','MN4'])
  def _dump_TidalRun(self):
    self.TidalRun.dump('./')


if __name__ == "__main__":
  _90DayTidal()
