import unittest
from AdcircPy import AdcircPy
from AdcircPyTests import AdcircPyEnvironment

class FrontEndTests(AdcircPyEnvironment, unittest.TestCase):
  def setUp(self):
    super(FrontEndTests, self).__init__()
    self.read_environment_variables()
  
  def test_read_mesh_coldstart(self):
    AdcircPy.read_mesh(fort14=self._os.getenv("FORT14_PATH"),
                       fort13=self._os.getenv("FORT13_PATH"),
                       fort15=self._os.getenv("FORT15_COLDSTART_PATH"))
    
  def test_read_mesh_hotstart(self): 
    AdcircPy.read_mesh(fort14=self._os.getenv("FORT14_PATH"),
                       fort13=self._os.getenv("FORT13_PATH"),
                       fort15=self._os.getenv("FORT15_HOTSTART_PATH"))

  def test_read_ascii_maxele(self):
    AdcircPy.read_output(self._os.getenv("MAXELE_ASCII_PATH"),
                         fort14=self._os.getenv("FORT14_PATH"))

  def test_read_ascii_maxele_no_fort14(self):
    with self.assertRaises(Exception) as context:
      AdcircPy.read_output(self._os.getenv("MAXELE_ASCII_PATH"))
    self.assertTrue('a fort.14 is required' in str(context.exception))

  def test_read_netcdf_maxele(self):
    AdcircPy.read_output(self._os.getenv("MAXELE_ASCII_PATH"),
                         fort14=self._os.getenv("FORT14_PATH"))

  def test_read_netcdf_maxele_no_fort14(self):
    """
    NOTE: This function is returning True
    but in reality the fort.14 file is being ignored.
    Need to collect boundary data from fort.14 before returning object.
    """
    AdcircPy.read_output(self._os.getenv("MAXELE_NC_PATH"),
                         fort14=self._os.getenv("FORT14_PATH"))



if __name__ == '__main__':
    unittest.main()

