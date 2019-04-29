# local imports
from AdcircPy.Mesh.Boundaries._BaseBoundary import _BaseBoundary

# unittest imports
import unittest

class CulvertBoundaries(_BaseBoundary):

    def __init__(self, *boundaries):
        super(CulvertBoundaries, self).__init__(*boundaries)

    def add_boundary(self, front_face, back_face, ):
        pass


class CulvertBoundariesTestCase(unittest.TestCase):

    def setUp(self):
        self.CulvertBoundaries = CulvertBoundaries

    def test_empty(self):
        self.CulvertBoundaries()