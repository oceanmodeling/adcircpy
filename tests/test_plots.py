import unittest
import matplotlib
if os.getenv('CIRCLECI') == 'true':
  matplotlib.use('Agg')
import matplotlib.pyplot as plt 
from AdcircPy import AdcircPy
from AdcircPyTests import AdcircPyEnvironment


class FrontEndTests(AdcircPyEnvironment, unittest.TestCase):
  def setUp(self):
    super(FrontEndTests, self).__init__()
    self.read_environment_variables()
    self.AdcircMesh = AdcircPy.read_mesh(fort14=self._os.getenv("FORT14_PATH"))

  def test_make_plot(self):
    self.AdcircMesh.make_plot()
    plt.show()
    plt.close(plt.gcf())


  def test_plot_trimesh(self):
    self.AdcircMesh.plot_trimesh()
    plt.show()
    plt.close(plt.gcf())


if __name__ == '__main__':
    unittest.main()

