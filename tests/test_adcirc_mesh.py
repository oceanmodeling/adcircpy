#! /usr/bin/env python
from pathlib import Path
import tempfile

import pytest as pytest

from adcircpy import AdcircMesh


@pytest.fixture
def nodes() -> {str: ((float, float), float)}:
    return {
        '1': ((0.0, 0.0), -5.0),
        '2': ((0.5, 0.0), -4.0),
        '3': ((1.0, 0.0), -3.0),
        '4': ((1.0, 1.0), -2.0),
        '5': ((0.0, 1.0), -1.0),
        '6': ((0.5, 1.5), 0.0),
        '7': ((0.33, 0.33), 1.0),
        '8': ((0.66, 0.33), 2.0),
        '9': ((0.5, 0.66), 3.0),
        '10': ((-1.0, 1.0), 4.0),
        '11': ((-1.0, 0.0), 5.0),
    }


@pytest.fixture
def elements() -> {str: [str]}:
    return {
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


@pytest.fixture
def boundaries() -> {int: {int: {str: [str]}}}:
    return {
        # "open" boundaries
        None: {0: {'indexes': ['10', '11', '1', '2']}, 1: {'indexes': ['2', '3', '4']},},
        # "land" boundaries
        0: {0: {'indexes': ['4', '6']}, 1: {'indexes': ['6', '5', '10']},},
        # "interior" boundary
        1: {0: {'indexes': ['7', '8', '9', '7']}},
    }


def test_triangles_only(nodes, elements):
    mesh = AdcircMesh(nodes, {id: geom for geom in elements.values() if len(geom) == 3},)

    assert isinstance(mesh, AdcircMesh)


def test_quads_only(nodes, elements):
    mesh = AdcircMesh(nodes, {id: geom for geom in elements.values() if len(geom) == 4},)

    assert isinstance(mesh, AdcircMesh)


def test_hybrid(nodes, elements):
    mesh = AdcircMesh(nodes, elements)

    assert isinstance(mesh, AdcircMesh)


def test_open(nodes, elements):
    temporary_file_handle = tempfile.NamedTemporaryFile()
    with open(temporary_file_handle.name, 'w') as temporary_file:
        temporary_file.write(f'\n{len(elements):d} {len(nodes):d}\n')
        for id, ((x, y), z) in nodes.items():
            temporary_file.write(f'{id} {x} {y} {z}\n')
        for id, geometry in elements.items():
            temporary_file.write(f'{id} {len(geometry)} {" ".join(idx for idx in geometry)}\n')

    mesh = AdcircMesh.open(temporary_file_handle.name)

    assert isinstance(mesh, AdcircMesh)


def test_make_plot(nodes, elements, mocker):
    mocker.patch('matplotlib.pyplot.show')

    mesh = AdcircMesh(nodes, elements)
    mesh.make_plot(
        show=True, extent=[0, 1, 0, 1], title='test', cbar_label='elevation [m]', vmax=0.0,
    )

    assert isinstance(mesh, AdcircMesh)


def test_make_plot_wet_only(mocker):
    mocker.patch('matplotlib.pyplot.show')

    nodes = {
        0: ((0.0, 0.0), 0.0),
        1: ((1.0, 0.0), -1.0),
        2: ((1.0, 1.0), -2.0),
        3: ((0.0, 1.0), -3.0),
        4: ((0.5, 1.5), -4.0),
    }
    elements = {
        0: [2, 4, 3],
        1: [0, 1, 2, 3],
    }

    mesh = AdcircMesh(nodes, elements)
    mesh.make_plot()

    assert isinstance(mesh, AdcircMesh)


def test_write(nodes, elements):
    mesh = AdcircMesh(nodes, elements)
    temporary_directory = tempfile.TemporaryDirectory()
    temporary_directory_path = Path(temporary_directory.name)

    mesh.write(Path(temporary_directory.name) / 'test_AdcircMesh.gr3',)
    mesh.write(
        Path(temporary_directory.name) / 'test_AdcircMesh.2dm', format='2dm',
    )

    with pytest.raises(Exception):
        mesh.write(
            temporary_directory_path / 'test_AdcircMesh.2dm', format='2dm',
        )

    with pytest.raises(ValueError):
        mesh.write(
            temporary_directory_path / 'test_AdcircMesh.txt', format='txt',
        )


def test_triplot(nodes, elements, boundaries, mocker):
    mocker.patch('matplotlib.pyplot.show')

    mesh = AdcircMesh(nodes, elements, boundaries)
    mesh.triplot()


def test_make_plot_flat_domain(nodes, elements, boundaries, mocker):
    mocker.patch('matplotlib.pyplot.show')

    nodes = {id: (coord, 0.0) for id, (coord, _) in nodes.items()}
    mesh = AdcircMesh(nodes, elements, boundaries)
    mesh.make_plot()
