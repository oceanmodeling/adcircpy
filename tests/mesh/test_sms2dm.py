#! /usr/bin/env python
import pathlib
import tempfile
import unittest

from adcircpy.mesh.sms2dm import reader, writer


class Sms2dmTestCase(unittest.TestCase):
    def setUp(self):
        self.nodes = {
            '1': ([0.0, 0.0], -99999.0),
            '2': ([0.5, 0.0], -99999.0),
            '3': ([1.0, 0.0], -99999.0),
            '4': ([1.0, 1.0], -99999.0),
            '5': ([0.0, 1.0], -99999.0),
            '6': ([0.5, 1.5], -99999.0),
            '7': ([0.33, 0.33], -99999.0),
            '8': ([0.66, 0.33], -99999.0),
            '9': ([0.5, 0.66], -99999.0),
            '10': ([-1.0, 1.0], -99999.0),
            '11': ([-1.0, 0.0], -99999.0),
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
            '10': ['5', '1', '7'],
        }

        self.boundaries = dict()

        self.boundaries[None] = {  # "open" boundaries
            0: {'indexes': ['10', '11', '1', '2']},
            1: {'indexes': ['2', '3', '4']},
        }

        self.boundaries[0] = {  # "land" boundaries
            0: {'indexes': ['4', '6']},
            1: {'indexes': ['6', '5', '10']},
        }

        self.boundaries[1] = {
            0: {'indexes': ['7', '8', '9', '7']}}  # "interior" boundary

    def test_write_read(self):
        grd = {
            'ND': self.nodes,
            'E3T': {id: indexes for id, indexes in self.elements.items() if
                    len(indexes) == 3},
            'E4Q': {id: indexes for id, indexes in self.elements.items() if
                    len(indexes) == 4},
            'boundaries': self.boundaries,
        }
        tmpdir = tempfile.TemporaryDirectory()
        tmpfile = pathlib.Path(tmpdir.name) / 'hgrid.2dm'
        writer(grd, pathlib.Path(tmpfile))
        grd.pop('boundaries')
        self.assertDictEqual(reader(pathlib.Path(tmpfile)), grd)

    def test_write_raise_exists(self):
        grd = {
            'ND': self.nodes,
            'E3T': {id: indexes for id, indexes in self.elements.items() if
                    len(indexes) == 3},
            'E4Q': {id: indexes for id, indexes in self.elements.items() if
                    len(indexes) == 4},
        }
        tmpdir = tempfile.TemporaryDirectory()
        tmpfile = pathlib.Path(tmpdir.name) / 'hgrid.2dm'
        writer(grd, pathlib.Path(tmpfile))
        self.assertRaises(Exception, writer, grd, pathlib.Path(tmpfile))


if __name__ == '__main__':
    unittest.main()
