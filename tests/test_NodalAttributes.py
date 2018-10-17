#! /usr/bin/env python
import unittest
from AdcircPy.Model import NodalAttributes

class TestNodalAttributes(unittest.TestCase):
  """
  unittest class for NodalAttributes class
  """
  def test_empty(self):
    NodalAttributes()

  def test_raise_wrong_attribute(self):
    with self.assertRaises(Exception) as context:
      NodalAttributes(spinup_attributes='None')
    self.assertTrue('not found in fort.13' in str(context.exception))
    
    

if __name__ == '__main__':
  unittest.main()