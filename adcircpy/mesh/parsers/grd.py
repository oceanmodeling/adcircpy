from collections import defaultdict
import numbers
import os
import pathlib
from typing import Dict, Iterable, TextIO, Union
import warnings

import numpy as np  # type: ignore[import]
from pyproj import CRS  # type: ignore[import]
from pyproj.exceptions import CRSError  # type: ignore[import]


def buffer_to_dict(buf: TextIO):
    description = buf.readline().strip()
    NE, NP = map(int, buf.readline().split())
    nodes = {}
    for _ in range(NP):
        line = buf.readline().strip('\n').split()
        # Gr3/fort.14 format cannot distinguish between a 2D mesh with one
        # vector value (e.g. velocity, which uses 2 columns) or a 3D mesh with
        # one scalar value. This is a design problem of the mesh format, which
        # renders it ambiguous, and the main reason why the use of fort.14/grd
        # formats is discouraged, in favor of UGRID.
        # Here, we assume the input mesh is strictly a 2D mesh, and the data
        # that follows is an array of values.
        if len(line[3:]) == 1:
            nodes[line[0]] = [(float(line[1]), float(line[2])), float(line[3])]
        else:
            nodes[line[0]] = [
                (float(line[1]), float(line[2])),
                [float(line[i]) for i in range(3, len(line[3:]))],
            ]
    elements = {}
    for _ in range(NE):
        line = buf.readline().split()
        elements[line[0]] = line[2:]
    # Assume EOF if NOPE is empty.
    try:
        NOPE = int(buf.readline().split()[0])
    except IndexError:
        return {'description': description, 'nodes': nodes, 'elements': elements}
    # let NOPE=-1 mean an ellipsoidal-mesh
    # reassigning NOPE to 0 until further implementation is applied.
    boundaries: Dict = defaultdict(dict)
    _bnd_id = 0
    buf.readline()
    while _bnd_id < NOPE:
        NETA = int(buf.readline().split()[0])
        _cnt = 0
        boundaries[None][_bnd_id] = dict()
        boundaries[None][_bnd_id]['node_id'] = list()
        while _cnt < NETA:
            boundaries[None][_bnd_id]['node_id'].append(buf.readline().split()[0].strip())
            _cnt += 1
        _bnd_id += 1
    NBOU = int(buf.readline().split()[0])
    _nbnd_cnt = 0
    buf.readline()
    while _nbnd_cnt < NBOU:
        npts, ibtype = buf.readline().split()[:2]
        _pnt_cnt = 0
        if ibtype not in boundaries:
            _bnd_id = 0
        else:
            _bnd_id = len(boundaries[ibtype])
        boundaries[ibtype][_bnd_id] = dict()
        boundaries[ibtype][_bnd_id]['node_id'] = list()
        if ibtype.endswith('3'):  # outflow
            boundaries[ibtype][_bnd_id]['barrier_height'] = list()
            boundaries[ibtype][_bnd_id]['supercritical_flow_coefficient'] = list()
        elif ibtype.endswith('4'):  # weir
            boundaries[ibtype][_bnd_id]['barrier_height'] = list()
            boundaries[ibtype][_bnd_id]['subcritical_flow_coefficient'] = list()
            boundaries[ibtype][_bnd_id]['supercritical_flow_coefficient'] = list()
        elif ibtype.endswith('5'):  # culvert
            boundaries[ibtype][_bnd_id]['barrier_height'] = list()
            boundaries[ibtype][_bnd_id]['subcritical_flow_coefficient'] = list()
            boundaries[ibtype][_bnd_id]['supercritical_flow_coefficient'] = list()
            boundaries[ibtype][_bnd_id]['cross_barrier_pipe_height'] = list()
            boundaries[ibtype][_bnd_id]['friction_factor'] = list()
            boundaries[ibtype][_bnd_id]['pipe_diameter'] = list()

        while _pnt_cnt < int(npts):
            line = buf.readline().split()
            if ibtype.endswith('3'):
                boundaries[ibtype][_bnd_id]['node_id'].append((line[0],))
                boundaries[ibtype][_bnd_id]['barrier_height'].append(float(line[1]))
                boundaries[ibtype][_bnd_id]['supercritical_flow_coefficient'].append(
                    float(line[2])
                )
            elif ibtype.endswith('4'):
                boundaries[ibtype][_bnd_id]['node_id'].append((line[0], line[1]))
                boundaries[ibtype][_bnd_id]['barrier_height'].append(float(line[2]))
                boundaries[ibtype][_bnd_id]['subcritical_flow_coefficient'].append(
                    float(line[3])
                )
                boundaries[ibtype][_bnd_id]['supercritical_flow_coefficient'].append(
                    float(line[4])
                )
            elif ibtype.endswith('5'):
                boundaries[ibtype][_bnd_id]['node_id'].append((line[0], line[1]))
                boundaries[ibtype][_bnd_id]['barrier_height'].append(float(line[2]))
                boundaries[ibtype][_bnd_id]['subcritical_flow_coefficient'].append(
                    float(line[3])
                )
                boundaries[ibtype][_bnd_id]['supercritical_flow_coefficient'].append(
                    float(line[4])
                )
                boundaries[ibtype][_bnd_id]['supercritical_flow_coefficient'].append(
                    float(line[4])
                )
                boundaries[ibtype][_bnd_id]['supercritical_flow_coefficient'].append(
                    float(line[4])
                )
                boundaries[ibtype][_bnd_id]['supercritical_flow_coefficient'].append(
                    float(line[4])
                )

            else:
                boundaries[ibtype][_bnd_id]['node_id'].append(line[0])
            _pnt_cnt += 1
        _nbnd_cnt += 1
    return {
        'description': description,
        'nodes': nodes,
        'elements': elements,
        'boundaries': boundaries,
    }


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
    for id, (coords, values) in nodes.items():
        if isinstance(values, numbers.Number):
            values = [values]
        line = [f'{id}']
        line.extend([f'{x:<.16E}' for x in coords])
        line.extend([f'{x:<.16E}' for x in values])
        out.append(' '.join(line))

    for id, element in elements.items():
        line = [f'{id}']
        line.append(f'{len(element)}')
        line.extend([f'{e}' for e in element])
        out.append(' '.join(line))

    # ocean boundaries
    if boundaries is not None:
        ocean_boundaries = boundaries.get(None, {})
        out.append(f'{len(ocean_boundaries):d} ' '! total number of ocean boundaries')
        # count total number of ocean boundaries
        _sum = 0
        for bnd in ocean_boundaries.values():
            _sum += len(bnd['node_id'])
        out.append(f'{int(_sum):d} ! total number of ocean boundary nodes')
        # write ocean boundary indexes
        for i, boundary in ocean_boundaries.items():
            out.append(
                f"{len(boundary['node_id']):d}" f' ! number of nodes for ocean_boundary_{i}'
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
            for bnd in boundaries[ibtype].values():
                _cnt += np.asarray(bnd['node_id']).shape[0]
    out.append(f'{_cnt:d} ! Total number of non-ocean boundary nodes')
    # all additional boundaries
    for ibtype, boundaries in boundaries.items():
        if ibtype is None:
            continue
        for id, boundary in boundaries.items():
            line = [
                f"{len(boundary['node_id']):d}",
                f'{ibtype}',
                f'! boundary {ibtype}:{id}',
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
    resource = pathlib.Path(resource)
    with open(resource, 'r') as stream:
        grd = buffer_to_dict(stream)
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
            f'File {str(resource)} does not contain CRS ' 'information and no CRS was given.'
        )
    if crs is not False:
        grd.update({'crs': crs})
    return grd


def write(grd, path, overwrite=False):
    path = pathlib.Path(path)
    if path.is_file() and not overwrite:
        raise FileExistsError('File exists, pass overwrite=True to allow overwrite.')
    with open(path, 'w') as f:
        f.write(to_string(**grd))
