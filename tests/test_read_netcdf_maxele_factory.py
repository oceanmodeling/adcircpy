#! /usr/bin/env python
import os
import sys
from AdcircPy import AdcircPy

class testReadOutputs(object):
  def __init__(self):
    self.read_environment_variables()
    self.test_read_maxele_ascii()

  def test_read_maxele_ascii(self):
    maxele = AdcircPy.read_output(self._os.getenv("MAXELE_NC_PATH"),
                         # fort14=self._os.getenv("FORT14_PATH")
                         )
    # maxele.make_plot()


  def read_environment_variables(self):
    with open(os.path.expanduser("~")+'/.adcpyrc') as lines:
      for line in lines:
        if 'export' in line and '#' not in line:
          line = line.strip('export \n').split('=')
          os.environ[line[0]] = line[1]
    self._os=os

if __name__ == '__main__':
  testReadOutputs()

