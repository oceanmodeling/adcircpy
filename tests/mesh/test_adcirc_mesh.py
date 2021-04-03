#! /usr/bin/env python
import pathlib
import tempfile
import unittest
from unittest.mock import patch


from adcircpy import AdcircMesh


class AdcircMeshTestCase(unittest.TestCase):

    def setUp(self):
        self.nodes = {
            '1': ((0., 0.), -5.),
            '2': ((.5, 0.), -4.),
            '3': ((1., 0.), -3.),
            '4': ((1., 1.), -2.),
            '5': ((0., 1.), -1.),
            '6': ((.5, 1.5), 0.),
            '7': ((.33, .33), 1.),
            '8': ((.66, .33), 2.),
            '9': ((.5, .66), 3.),
            '10': ((-1., 1.), 4.),
            '11': ((-1., 0.), 5.),
        }
        self.elements = {
            '1': ['5', '7', '9'],
            '2': ['1', '2', '7'],
            '3': ['2', '3', '8'],
            '4': ['8', '7', '2'],
            '5': ['3', '4', '8'],
            '6': ['4', '9', '8'],
            '7': ['4', '6', '5'],
            '8': ['5', '10', '11', '1'],
            '9': ['9', '4', '5'],
            '10': ['5', '1', '7']
        }

        self.boundaries = dict()

        self.boundaries[None] = {  # "open" boundaries
            0: {'indexes': ['10', '11', '1', '2']},
            1: {'indexes': ['2', '3', '4']}
        }

        self.boundaries[0] = {  # "land" boundaries
            0: {'indexes': ['4', '6']},
            1: {'indexes': ['6', '5', '10']}
        }

        self.boundaries[1] = {  # "interior" boundary
            0: {'indexes': ['7', '8', '9', '7']}
        }

        self.grd = {
            'nodes': self.nodes,
            'elements': self.elements,
            'boundaries': self.boundaries,
            'description': 'gr3_unittest'
        }

    def test_triangles_only(self):
        self.assertIsInstance(
            AdcircMesh(
                self.nodes,
                {id: geom for geom in self.elements.values() if len(geom) == 3}
            ),
            AdcircMesh
        )

    def test_quads_only(self):
        self.assertIsInstance(
            AdcircMesh(
                self.nodes,
                {id: geom for geom in self.elements.values() if len(geom) == 4}
            ),
            AdcircMesh
        )

    def test_hybrid(self):
        self.assertIsInstance(AdcircMesh(self.nodes, self.elements),
                              AdcircMesh)

    def test_open(self):
        tmpfile = tempfile.NamedTemporaryFile()
        with open(tmpfile.name, 'w') as f:
            f.write(f'\n{len(self.elements):d} {len(self.nodes):d}\n')
            for id, ((x, y), z) in self.nodes.items():
                f.write(f'{id} {x} {y} {z}\n')
            for id, geom in self.elements.items():
                f.write(f'{id} {len(geom)} {" ".join(idx for idx in geom)}\n')
        self.assertIsInstance(AdcircMesh.open(tmpfile.name), AdcircMesh)

    @patch('matplotlib.pyplot.show')
    def test_make_plot(self, mock):
        h = AdcircMesh(self.nodes, self.elements)
        h.make_plot(
            show=True,
            extent=[0, 1, 0, 1],
            title='test',
            cbar_label='elevation [m]',
            vmax=0.
        )
        self.assertIsInstance(h, AdcircMesh)

    @patch('matplotlib.pyplot.show')
    def test_make_plot_wet_only(self, mock):
        nodes = {
            0: ((0., 0.), 0.),
            1: ((1., 0.), -1.),
            2: ((1., 1.), -2.),
            3: ((0., 1.), -3.),
            4: ((0.5, 1.5), -4.),
        }
        elements = {
            0: [2, 4, 3],
            1: [0, 1, 2, 3],
        }
        h = AdcircMesh(nodes, elements)
        h.make_plot()
        self.assertIsInstance(h, AdcircMesh)

    def test_write(self):
        h = AdcircMesh(self.nodes, self.elements)
        tmpdir = tempfile.TemporaryDirectory()
        h.write(pathlib.Path(tmpdir.name) / 'test_AdcircMesh.gr3')
        h.write(pathlib.Path(tmpdir.name) / 'test_AdcircMesh.2dm', format='2dm')
        self.assertRaises(Exception, h.write,
                          pathlib.Path(tmpdir.name) / 'test_AdcircMesh.2dm',
                          format='2dm')
        self.assertRaises(ValueError, h.write,
                          pathlib.Path(tmpdir.name) / 'test_AdcircMesh.txt',
                          format='txt')

    def test_triplot(self):
        h = AdcircMesh(self.nodes, self.elements, self.boundaries)
        h.triplot()

    @patch('matplotlib.pyplot.show')
    def test_make_plot_flat_domain(self, mock):
        nodes = {id: (coord, 0.) for id, (coord, _) in self.nodes.items()}
        h = AdcircMesh(nodes, self.elements, self.boundaries)
        h.make_plot()


if __name__ == '__main__':
    unittest.main()
