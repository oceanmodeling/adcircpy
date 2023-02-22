#! /usr/bin/env python

import pytest

from adcircpy.mesh import AdcircMesh

# noinspection PyUnresolvedReferences
from tests import shinnecock_mesh_directory


def test_mesh_make_plot(shinnecock_mesh_directory):
    mesh = AdcircMesh.open(shinnecock_mesh_directory / 'fort.14', crs=4326)
    mesh.make_plot()


def test_mesh_get_land_boundaries(shinnecock_mesh_directory):
    mesh = AdcircMesh.open(shinnecock_mesh_directory / 'fort.14', crs=4326)
    mesh.land_boundaries.gdf


def test_mesh_get_ocean_boundaries(shinnecock_mesh_directory):
    mesh = AdcircMesh.open(shinnecock_mesh_directory / 'fort.14', crs=4326)
    mesh.ocean_boundaries.gdf
