#! /usr/bin/env python
import os
import sys
sys.path.append('..')
from AdcircPy import AdcircPy

class testFort15(object):
  def __init__(self):
    self._getenv()
    self._test_coldstart()
    self._test_hotstart()

  def _getenv(self):
    with open(os.path.expanduser("~")+'/.adcpyrc') as lines:
      for line in lines:
        if 'export' in line and '#' not in line:
          line = line.strip('export \n').split('=')
          os.environ[line[0]] = line[1]

  def _test_coldstart(self):
    AdcircPy.read_mesh(fort14=os.getenv("ADCIRC_MESH_PATH"),
                       fort15=os.getenv("ADCIRC_FORT15_COLDSTART_PATH"))

  def _test_hotstart(self):
    AdcircPy.read_mesh(fort14=os.getenv("ADCIRC_MESH_PATH"),
                       fort15=os.getenv("ADCIRC_FORT15_HOTSTART_PATH"))

if __name__ == '__main__':
  testFort15()

