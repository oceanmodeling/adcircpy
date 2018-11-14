#! /usr/bin/env python
from AdcircPyTests import AdcircPyEnvironment
import unittest
import matplotlib
import matplotlib.pyplot as plt 
from AdcircPy import AdcircPy

class FrontEndTests(AdcircPyEnvironment, unittest.TestCase):
  def setUp(self):
    super(FrontEndTests, self).__init__()
    self.read_environment_variables()
 
  def test_read_harmonic_constituents_ascii(self):
    AdcircPy.read_output(self.os.getenv('HARM_CONST_ASCII'), fort15=self.os.getenv('HARM_CONST_ASCII_FORT15'))

  def test_read_harmonic_constituents_ascii_no_fort15(self):
    with self.assertRaises(Exception) as context:
      AdcircPy.read_output(self.os.getenv('HARM_CONST_ASCII'))
    self.assertTrue('fort.15 needs to be provided for fort.51' in str(context.exception))

  def test_read_harmonic_constituents_netcdf(self):
    AdcircPy.read_output(self.os.getenv('HARM_CONST_NETCDF'))

if __name__ == '__main__':
    unittest.main()

