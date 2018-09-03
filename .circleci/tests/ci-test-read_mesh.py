import unittest
from AdcircPy import AdcircPy
class testFrontEnd(AdcircPy, unittest.TestCase):
    def test_read_mesh(self):
        AdcircPy.read_mesh(fort14='fort.14',
                           fort13='fort.13',
                           fort15='fort.15')
if __name__ == '__main__':
    unittest.main()

