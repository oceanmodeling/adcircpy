import unittest
from datetime import datetime, timedelta
from AdcircPy import AdcircPy
from AdcircPyTests import AdcircPyEnvironment


class GenerateTidalRunTest(AdcircPyEnvironment, unittest.TestCase):
  def setUp(self):
    super(GenerateTidalRunTest, self).__init__()
    self.read_environment_variables()
  
  def test_generate_tidal_run_from_AdcircMesh(self):
    mesh = AdcircPy.read_mesh(fort14=self._os.getenv("FORT14_PATH"),
                              fort13=self._os.getenv("FORT13_PATH"))
    spinup_date = datetime(2013, 8, 1, 0, 0, 0)
    start_time = spinup_date + timedelta(days=15)
    end_time   = spinup_date + timedelta(days=90)
    tidal_run  = mesh.generate_tidal_run(start_time, end_time,
                                        # constituents=constituents
                                        spinup_date=spinup_date
                                        )
    
  

if __name__ == '__main__':
    unittest.main()