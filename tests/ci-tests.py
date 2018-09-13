import unittest
from AdcircPy import AdcircPy

class tests(AdcircPy, unittest.TestCase):
  
  def test_read_mesh_coldstart(self):
    AdcircPy.read_mesh(fort14='fort.14',
                       fort13='fort.13',
                       fort15='fort.15.coldstart')
    
  def test_read_mesh_hotstart(self): 
    AdcircPy.read_mesh(fort14='fort.14',
                       fort13='fort.13',
                       fort15='fort.15.hotstart')

  def test_read_ascii_maxele(self):
    AdcircPy.read_output('maxele.63', fort14='fort.14')

  def test_read_ascii_maxele_no_fort14(self):
    with self.assertRaises(Exception) as context:
      AdcircPy.read_output('maxele.63')
    self.assertTrue('a fort.14 is required' in str(context.exception))


if __name__ == '__main__':
    unittest.main()

