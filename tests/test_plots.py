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
    self.AdcircMesh = AdcircPy.read_mesh(fort14=self.os.getenv("FORT14"))

  def test_make_plot(self):
    self.AdcircMesh.make_plot(Trimesh=True)
    plt.show()
    plt.close(plt.gcf())


if __name__ == '__main__':
    unittest.main()

