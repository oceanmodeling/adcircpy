#! /usr/bin/env python

import pytest

from adcircpy import AdcircMesh
from tests import OUTPUT_DIRECTORY


@pytest.fixture
def nodes() -> {int: ((float, float), float)}:
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
def elements() -> {int: [int]}:
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
def boundaries() -> {int: {int: {str: [int]}}}:
    return {
        # "open" boundaries
        None: {0: {'indexes': ['10', '11', '1', '2']}, 1: {'indexes': ['2', '3', '4']}},
        # "land" boundaries
        0: {0: {'indexes': ['4', '6']}, 1: {'indexes': ['6', '5', '10']}},
        # "interior" boundary
        1: {0: {'indexes': ['7', '8', '9', '7']}},
    }


@pytest.fixture
def fort14(elements, nodes) -> str:
    lines = [
        f'\n{len(elements):d} {len(nodes):d}',
        *(f'{id} {x} {y} {z}' for id, ((x, y), z) in nodes.items()),
        *(
            f'{id} {len(geometry)} {" ".join(idx for idx in geometry)}'
            for id, geometry in elements.items()
        ),
    ]
    return '\n'.join(lines)


@pytest.fixture
def wet_nodes() -> {int: ((float, float), float)}:
    return {
        0: ((0.0, 0.0), 0.0),
        1: ((1.0, 0.0), -1.0),
        2: ((1.0, 1.0), -2.0),
        3: ((0.0, 1.0), -3.0),
        4: ((0.5, 1.5), -4.0),
    }


@pytest.fixture
def wet_elements() -> {int: [int]}:
    return {
        0: [2, 4, 3],
        1: [0, 1, 2, 3],
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


def test_open(fort14):
    output_directory = OUTPUT_DIRECTORY / 'test_open'

    if not output_directory.exists():
        output_directory.mkdir(parents=True, exist_ok=True)

    mesh_filename = output_directory / 'fort.14'

    with open(mesh_filename, 'w') as temporary_file:
        temporary_file.write(fort14)

    mesh = AdcircMesh.open(mesh_filename)

    assert isinstance(mesh, AdcircMesh)


def test_make_plot(nodes, elements, mocker):
    mocker.patch('matplotlib.pyplot.show')

    mesh = AdcircMesh(nodes, elements)
    mesh.make_plot(
        show=True, extent=[0, 1, 0, 1], title='test', cbar_label='elevation [m]', vmax=0.0,
    )

    assert isinstance(mesh, AdcircMesh)


def test_make_plot_wet_only(wet_nodes, wet_elements, mocker):
    mesh = AdcircMesh(wet_nodes, wet_elements)

    mocker.patch('matplotlib.pyplot.show')
    mesh.make_plot()

    assert isinstance(mesh, AdcircMesh)


def test_write(nodes, elements):
    output_directory = OUTPUT_DIRECTORY / 'test_write'

    if not output_directory.exists():
        output_directory.mkdir(parents=True, exist_ok=True)

    mesh = AdcircMesh(nodes, elements)

    mesh.write(output_directory / 'test_AdcircMesh.gr3', overwrite=True)
    mesh.write(output_directory / 'test_AdcircMesh.2dm', format='2dm', overwrite=True)

    with pytest.raises(FileExistsError):
        mesh.write(output_directory / 'test_AdcircMesh.2dm', format='2dm')

    with pytest.raises(ValueError):
        mesh.write(output_directory / 'test_AdcircMesh.txt', format='txt', overwrite=True)


def test_triplot(nodes, elements, boundaries, mocker):
    mesh = AdcircMesh(nodes, elements, boundaries)

    mocker.patch('matplotlib.pyplot.show')
    mesh.triplot()


def test_make_plot_flat_domain(nodes, elements, boundaries, mocker):
    nodes = {id: (coord, 0.0) for id, (coord, _) in nodes.items()}
    mesh = AdcircMesh(nodes, elements, boundaries)

    mocker.patch('matplotlib.pyplot.show')
    mesh.make_plot()
