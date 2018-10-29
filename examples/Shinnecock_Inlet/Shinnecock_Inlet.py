#! /usr/bin/env python
from datetime import datetime, timedelta
from AdcircPy import ElevationStationsOutput as ESO
from AdcircPy import AdcircPy


class ShinnecockInletExample(object):
  def __init__(self):
    self._read_mesh()
    self._init_dates()
    self._init_forcing_constituents()
    self._init_ElevationStationsOutput()
    self._generate_TidalRun()
    self._write_files()

  def _read_mesh(self):
    self.ShinnecockInletMesh = AdcircPy.read_mesh('fort.14')

  def _init_dates(self):
    # Select start and end date, and optionally a spinup date.
    self.spinup_date = datetime(2015, 12, 14)
    self.start_date = self.spinup_date + timedelta(days=2)
    self.end_date = self.start_date + timedelta(days=5)

  def _init_forcing_constituents(self):
    # Select constituents to force.
    self.constituents = ['M2', 'N2', 'S2', 'K1', 'O1']

  def _init_ElevationStationsOutput(self):
    self.ESO = ESO.from_fort15('./fort.15.elevation_stations')

  def _generate_TidalRun(self):
    # Generate tidal run
    self.TidalRun = self.ShinnecockInletMesh.TidalRun(self.start_date, self.end_date,
                                                      spinup_date=self.spinup_date,
                                                      constituents=self.constituents,
                                                      DTDP=6.,
                                                      ElevationStationsOutput=self.ESO)

  def _write_files(self):
    # write pair of fort.15's to current directory.
    self.TidalRun.dump('./')

if __name__ == '__main__':
  ShinnecockInletExample()