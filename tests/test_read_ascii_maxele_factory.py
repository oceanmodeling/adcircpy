#! /usr/bin/env python
import unittest
import matplotlib
matplotlib.use('Agg')
from AdcircPyTests import AdcircPyEnvironment
from AdcircPy import AdcircPy

class testReadMaxeleAscii(AdcircPyEnvironment, unittest.TestCase):
  def setUp(self):
    self.read_environment_variables()

  def test_read_maxele_ascii(self):
    maxele = AdcircPy.read_output(self._os.getenv("MAXELE_ASCII_PATH"),
                                  fort14=self._os.getenv("FORT14_PATH"))
    maxele.make_plot()

class testReadOutputs(AdcircPyEnvironment, unittest.TestCase):
  def setUp(self):
    self.read_environment_variables()

  def test_read_maxele_ascii(self):
    maxele = AdcircPy.read_output(self._os.getenv("MAXELE_NC_PATH"))
    maxele.make_plot()

if __name__ == '__main__':
  unittest.main()

