#! /usr/bin/env python
"""
Example program for generating the official Sandy run for HSOFS.
The fort.14 and fort.13 files provided should correspond to HSOFS mesh.
The fort.15 file provided needs only to contain a list of elevation stations
preceeded by the keyword NOUTE anywhere in the file. These elevation stations
will be copied over to the newly generated fort.15. Elevation stations is the only
information parsed from the fort.15. Everything else (forcings, et al)
is generated anew.

fort.14 and fort.13 are large files and therefore they not provided in the example
directory. However, these files are publicly available and may be requested on the 
ADCIRC list-serv. The program should still work with any arbitrary mesh, so any other
mesh can in principle be provided for testing purposes.

The main purpose of this program is to show the users how a Best Track run can be
setup using this package, and should serve as an example so the users can build other
test cases.
"""

import os
from datetime import datetime, timedelta
from AdcircPy import AdcircPy
from AdcircPy import ElevationGlobalOutput as EGO
from AdcircPy import ElevationStationsOutput as ESO

MESH_PATH='./fort.14'
FORT13_PATH='./fort.13'
FORT15_PATH='./fort.15.original'

class SandyOfficial(object):

  def __init__(self):
    self._init_mesh()
    self._init_dates()
    self._init_elevation_station_output_request()
    self._init_elevation_global_output_request()
    self._init_BestTrackRun()
    self._dump_BestTrackRun()

  def _init_mesh(self):
    self.AdcircMesh = AdcircPy.read_mesh(fort14=MESH_PATH,
                                         fort13=FORT13_PATH)

  def _init_dates(self):
    self.spinup_date = datetime(2012, 10, 11, 0)
    self.start_time  = self.spinup_date + timedelta(days=15)
    self.end_time    = self.spinup_date + timedelta(days=19.25)
  
  def _init_elevation_station_output_request(self):
    self.ElevationStationsOutput = ESO.from_fort15(FORT15_PATH, sampling_frequency=timedelta(minutes=6))

  def _init_elevation_global_output_request(self):
    self.ElevationGlobalOutput = EGO(sampling_frequency=timedelta(minutes=15))

  def _init_BestTrackRun(self):
    self.BestTrackRun = self.AdcircMesh.BestTrackRun('AL182012', self.start_time, self.end_time,
                                        ElevationStationsOutput=self.ElevationStationsOutput,
                                        ElevationGlobalOutput=self.ElevationGlobalOutput,
                                        constituents=['K1', 'O1', 'P1', 'Q1',
                                                      'N2', 'M2', 'S2', 'K2',
                                                      'Mf', 'Mm', 'M4', 'MS4', 'MN4'])

  def _dump_BestTrackRun(self):
    self.BestTrackRun.dump('./official')

if __name__ == "__main__":
  SandyOfficial()
