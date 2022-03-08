from io import StringIO
import logging
import os
import pathlib
from typing import Iterable, Union
import warnings

import numpy as np  # type: ignore[import]
import pandas
from pyproj import CRS  # type: ignore[import]
from pyproj.exceptions import CRSError  # type: ignore[import]


def read_fort14(filename: os.PathLike):
    output = {'description': None, 'nodes': None, 'elements': None}
    with open(filename, 'r') as file:
        output['description'] = file.readline().strip()
        num_elements, num_nodes = [int(value) for value in file.readline().split()]

        # Gr3/fort.14 format cannot distinguish between a 2D mesh with one vector value (e.g. velocity, which uses 2 columns) or a 3D mesh with one scalar value. This is a design problem of the mesh format, which renders it ambiguous, and the main reason why the use of fort.14/grd formats is discouraged, in favor of UGRID. Here, we assume the input mesh is strictly a 2D mesh, and the data that follows is an array of values.
        with StringIO('\n'.join(file.readline() for _ in range(num_nodes))) as nodes_stream:
            output['nodes'] = pandas.read_csv(
                nodes_stream,
                delim_whitespace=True,
                index_col=0,
                names=[
                    'id',
                    'x',
                    'y',
                    *(f'value_{value_index}' for value_index in range(1, 5)),
                ],
            )

        with StringIO(
            '\n'.join(file.readline() for _ in range(num_elements))
        ) as elements_stream:
            output['elements'] = pandas.read_csv(
                elements_stream,
                delim_whitespace=True,
                index_col=0,
                names=['id', 'n', *(f'node_{node_index}' for node_index in range(1, 10))],
            )

        boundaries = {}

        # return on EOF
        try:
            num_open_boundaries = int(file.readline().split(maxsplit=1)[0])
        except IndexError:
            return output
        # num_open_boundary_nodes
        file.readline()

        # let NOPE=-1 mean an ellipsoidal-mesh
        # reassigning NOPE to 0 until further implementation is applied.
        if num_open_boundaries > 0:
            boundaries[None] = []
            for _ in range(num_open_boundaries):
                num_boundary_nodes = int(file.readline().split(maxsplit=1)[0])
                boundary = {
                    'node_id': [
                        int(file.readline().split(maxsplit=1)[0].strip())
                        for _ in range(num_boundary_nodes)
                    ]
                }
                boundaries[None].append(boundary)

        # return on EOF
        try:
            num_land_boundaries = int(file.readline().split(maxsplit=1)[0])
        except IndexError:
            return output
        # num_land_boundary_nodes
        file.readline()

        if num_land_boundaries > 0:
            for _ in range(num_land_boundaries):
                num_boundary_nodes, boundary_type = file.readline().split(maxsplit=2)[:2]
                num_boundary_nodes = int(num_boundary_nodes)
                boundary = {'node_id': []}

                last_digit = int(boundary_type[-1])

                if last_digit == 3:  # outflow
                    boundary.update(
                        {'barrier_height': [], 'supercritical_flow_coefficient': [],}
                    )
                elif last_digit == 4:  # weir
                    boundary.update(
                        {
                            'barrier_height': [],
                            'subcritical_flow_coefficient': [],
                            'supercritical_flow_coefficient': [],
                        }
                    )
                elif last_digit == 5:  # culvert
                    boundary.update(
                        {
                            'barrier_height': [],
                            'subcritical_flow_coefficient': [],
                            'supercritical_flow_coefficient': [],
                            'cross_barrier_pipe_height': [],
                            'friction_factor': [],
                            'pipe_diameter': [],
                        }
                    )

                for _ in range(num_boundary_nodes):
                    line = file.readline().split()
                    if last_digit == 3:
                        node = {
                            'node_id': (line[0],),
                            'barrier_height': float(line[1]),
                            'supercritical_flow_coefficient': float(line[2]),
                        }
                    elif last_digit == 4:
                        node = {
                            'node_id': (line[0], line[1]),
                            'barrier_height': float(line[2]),
                            'subcritical_flow_coefficient': float(line[3]),
                            'supercritical_flow_coefficient': float(line[4]),
                        }
                    elif last_digit == 5:
                        node = {
                            'node_id': (line[0], line[1]),
                            'barrier_height': float(line[2]),
                            'subcritical_flow_coefficient': float(line[3]),
                            'supercritical_flow_coefficient': float(line[4]),
                            'cross_barrier_pipe_height': float(line[5]),
                            'friction_factor': float(line[6]),
                            'pipe_diameter': float(line[7]),
                        }
                    else:
                        node = {
                            'node_id': line[0],
                        }

                    for key, value in node.items():
                        boundary[key].append(value)

                if boundary_type not in boundaries:
                    boundaries[boundary_type] = []

                boundaries[boundary_type].append(boundary)

        output['boundaries'] = boundaries
        return output


def to_string(description, nodes, elements, boundaries=None, crs=None):
    """
    must contain keys:
        description
        vertices
        elements
        vertex_id
        element_id
        values
        boundaries (optional)
            indexes
    """
    NE, NP = len(elements), len(nodes)
    out = [f'{description}', f'{NE} {NP}']

    # TODO: Probably faster if using np.array2string
    def float_format(value: float):
        return f'{value:<.16E}'

    out.append(
        nodes.to_string(
            None, header=False, index_names=False, justify='left', float_format=float_format
        )
    )
    out.append(
        elements.to_string(
            None,
            header=False,
            index_names=False,
            justify='left',
            float_format=float_format,
            na_rep='',
        )
    )

    # ocean boundaries
    if boundaries is not None:
        ocean_boundaries = boundaries.get(None, {})
        out.append(f'{len(ocean_boundaries):d} ' '! total number of ocean boundaries')
        # count total number of ocean boundaries
        _sum = 0
        for bnd in ocean_boundaries:
            _sum += len(bnd['node_id'])
        out.append(f'{int(_sum):d} ! total number of ocean boundary nodes')
        # write ocean boundary indexes
        for i, boundary in enumerate(ocean_boundaries):
            out.append(
                f'{len(boundary["node_id"]):d} ! number of nodes for ocean_boundary_{i + 1}'
            )
            for idx in boundary['node_id']:
                out.append(f'{idx}')
    else:
        out.append('0 ! total number of ocean boundaries')
        out.append('0 ! total number of ocean boundary nodes')
    # remaining boundaries
    boundaries = {} if boundaries is None else boundaries
    _cnt = 0
    for key in boundaries:
        if key is not None:
            for bnd in boundaries[key]:
                _cnt += 1
    out.append(f'{_cnt:d}  ! total number of non-ocean boundaries')
    # count remaining boundary nodes
    _cnt = 0
    for ibtype in boundaries:
        if ibtype is not None:
            for bnd in boundaries[ibtype]:
                _cnt += np.asarray(bnd['node_id']).shape[0]
    out.append(f'{_cnt:d} ! Total number of non-ocean boundary nodes')
    # all additional boundaries
    for ibtype, type_boundaries in boundaries.items():
        if ibtype is None:
            continue
        for index, boundary in enumerate(type_boundaries):
            line = [
                f"{len(boundary['node_id']):d}",
                f'{ibtype}',
                f'! boundary {ibtype}:{index + 1}',
            ]
            out.append(' '.join(line))
            for i, node_id in enumerate(boundary['node_id']):
                if isinstance(node_id, Iterable) and not isinstance(node_id, str):
                    line = [' '.join([x for x in list(node_id)])]
                else:
                    line = [node_id]
                if ibtype.endswith('3'):  # outflow
                    line.append(f'{boundary["barrier_height"][i]:.16e}')
                    line.append(f'{boundary["supercritical_flow_coefficient"][i]:.16e}')
                elif ibtype.endswith('4'):  # weir
                    line.append(f'{boundary["barrier_height"][i]:.16e}')
                    line.append(f'{boundary["subcritical_flow_coefficient"][i]:.16e}')
                    line.append(f'{boundary["supercritical_flow_coefficient"][i]:.16e}')
                elif ibtype.endswith('5'):  # culvert
                    line.append(f'{boundary["barrier_height"][i]:.16e}')
                    line.append(f'{boundary["subcritical_flow_coefficient"][i]:.16e}')
                    line.append(f'{boundary["supercritical_flow_coefficient"][i]:.16e}')
                    line.append(f'{boundary["cross_barrier_pipe_height"][i]:.16e}')
                    line.append(f'{boundary["friction_factor"][i]:.16e}')
                    line.append(f'{boundary["pipe_diameter"][i]:.16e}')
                out.append(' '.join(line))
    return '\n'.join(out)


def read(resource: Union[str, os.PathLike], boundaries: bool = True, crs=True):
    """Converts a file-like object representing a grd-formatted unstructured
    mesh into a python dictionary:

    Args:
        resource: Path to file on disk or file-like object such as
            :class:`io.StringIO`
    """

    if not isinstance(resource, pathlib.Path):
        resource = pathlib.Path(resource)
    grd = read_fort14(resource)
    if boundaries is False:
        grd.pop('boundaries', None)
    if crs is True:
        crs = None
    if crs is None:
        for try_crs in grd['description'].split():
            try:
                crs = CRS.from_user_input(try_crs)
                break
            except CRSError:
                pass
    if crs is None:
        warnings.warn(
            f'File does not contain CRS information and no CRS was given: "{str(resource)}"'
        )
    if crs is not False:
        grd.update({'crs': crs})
    return grd


def write(grd, path, overwrite=False):
    if not isinstance(path, pathlib.Path):
        path = pathlib.Path(path)
    if overwrite or not path.exists():
        with open(path, 'w') as f:
            f.write(to_string(**grd))
    else:
        logging.debug(f'skipping existing file "{path}"')
