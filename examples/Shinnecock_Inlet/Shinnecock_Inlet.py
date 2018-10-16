#! /usr/bin/env python
from AdcircPy import AdcircPy
from datetime import datetime, timedelta

# Read mesh
Mesh = AdcircPy.read_mesh('fort.14')

# Select start and end date, and optionally a spinup date.
spinup_date = datetime(2015, 12, 14)
start_date = spinup_date + timedelta(days=2)
end_date = start_date + timedelta(days=5)

# Select constituents to force.
constituents = ['M2', 'N2', 'S2', 'K1', 'O1']

# Generate tidal run
TidalRun = Mesh.TidalRun(start_date, end_date, spinup_date=spinup_date,
                           constituents=constituents, DTDP=6.)

# write files to current directory.
TidalRun.dump('./')