from collections import defaultdict
import pathlib

import numpy as np


def reader(path):
    grd = dict()
    grd['nodes'] = defaultdict(list)
    grd['elements'] = defaultdict(list)
    with open(pathlib.Path(path), 'r') as f:
        grd['description'] = f.readline().strip('\n')
        NE, NP = map(int, f.readline().split())
        for i in range(NP):
            id, x, y, z = f.readline().split()
            grd['nodes'][id] = ((float(x), float(y)), float(z))
        for i in range(NE):
            geom = f.readline().split()
            grd['elements'][geom[0]] = [x for x in geom[2:]]
        # Assume EOF if NOPE is empty.
        try:
            NOPE = int(f.readline().split()[0])
        except IndexError:
            return grd
        # let NOPE=-1 mean an ellipsoidal-mesh
        # reassigning NOPE to 0 until further implementation is applied.
        grd['boundaries'] = defaultdict(dict)
        _bnd_id = 0
        f.readline()  # Not used.
        while _bnd_id < NOPE:
            NETA = int(f.readline().split()[0])
            _cnt = 0
            grd['boundaries'][None][_bnd_id] = dict()
            grd['boundaries'][None][_bnd_id]['indexes'] = list()
            while _cnt < NETA:
                grd['boundaries'][None][_bnd_id]['indexes'].append(
                    f.readline().split()[0].strip())
                _cnt += 1
            _bnd_id += 1
        NBOU = int(f.readline().split()[0])
        _nbnd_cnt = 0
        f.readline()  # not used
        while _nbnd_cnt < NBOU:
            npts, ibtype = map(int, f.readline().split()[:2])
            _pnt_cnt = 0
            if ibtype not in grd['boundaries']:
                _bnd_id = 0
            else:
                _bnd_id = len(grd['boundaries'][ibtype])
            grd['boundaries'][ibtype][_bnd_id] = dict()
            grd['boundaries'][ibtype][_bnd_id]['indexes'] = list()
            while _pnt_cnt < npts:
                line = f.readline().split()
                grd['boundaries'][ibtype][_bnd_id]['indexes'].append(line[0])
                _pnt_cnt += 1
            _nbnd_cnt += 1
    # typecast defaultdict to regular dict
    for key, value in grd.items():
        if key != 'description':
            grd[key] = dict(value)
    return grd


def writer(grd, path, overwrite=False):
    path = pathlib.Path(path)
    if path.is_file() and not overwrite:
        msg = 'File exists, pass overwrite=True to allow overwrite.'
        raise Exception(msg)
    with open(path, 'w') as f:
        f.write(string(grd))
    return 0  # for unittests


def string(grd):
    f = graph(grd)
    if 'boundaries' in grd.keys():
        f += boundaries(grd)
    return f


def graph(grd):
    f = ""
    f += f"{grd['description']}\n"
    f += f"{len(grd['elements'])} "
    f += f"{len(grd['nodes'])}\n"
    # TODO: Make faster using np.array2string
    for id, ((x, y), z) in grd['nodes'].items():
        f += f"{id} "
        f += f"{x:<.16E} "
        f += f"{y:<.16E} "
        f += f"{z:<.16E}\n"
    for id, geom in grd['elements'].items():
        f += f"{id} "
        f += f"{len(geom):d} "
        for idx in geom:
            f += f"{idx} "
        f += "\n"
    return f


def boundaries(grd):
    f = ""
    # ocean boundaries
    if None in grd['boundaries']:
        f += f"{len(grd['boundaries'][None]):d} "
        f += "! total number of ocean boundaries\n"
        # count total number of ocean boundaries
        _sum = 0
        for bnd in grd['boundaries'][None].values():
            _sum += len(bnd['indexes'])
        f += f"{int(_sum):d} ! total number of ocean boundary nodes\n"
        # write ocean boundary indexes
        for i, boundary in grd['boundaries'][None].items():
            f += f"{len(boundary['indexes']):d}"
            f += f" ! number of nodes for ocean_boundary_{i}\n"
            for idx in boundary['indexes']:
                f += f"{idx}\n"
    else:
        f += "0 ! total number of ocean boundaries\n"
        f += "0 ! total number of ocean boundary nodes\n"
    # remaining boundaries
    _cnt = 0
    for key in grd['boundaries']:
        if key is not None:
            for bnd in grd['boundaries'][key]:
                _cnt += 1
    f += f"{_cnt:d}  ! total number of non-ocean boundaries\n"
    # count remaining boundary nodes
    _cnt = 0
    for ibtype in grd['boundaries']:
        if ibtype is not None:
            for bnd in grd['boundaries'][ibtype].values():
                _cnt += np.asarray(bnd['indexes']).size
    f += f"{_cnt:d} ! Total number of non-ocean boundary nodes\n"
    # all additional boundaries
    for ibtype, boundaries in grd['boundaries'].items():
        if ibtype is None:
            continue
        for id, boundary in boundaries.items():
            f += f"{len(boundary['indexes']):d} "
            f += f"{ibtype} "
            f += f"! boundary {ibtype}:{id}\n"
            for idx in boundary['indexes']:
                f += f"{idx}\n"
    return f


def euclidean_mesh(grd):
    # from geomesh.mesh.hgrid import Hgrid
    # cast gr3 inputs into a geomesh structure format
    coords = {id: (x, y) for id, ((x, y), value) in grd['nodes'].items()}
    triangles = {id: geom for id, geom in grd['elements'].items()
                 if len(geom) == 3}
    quads = {id: geom for id, geom in grd['elements'].items()
             if len(geom) == 4}
    values = [-value for coord, value in grd['nodes'].values()]
    description = grd['description']  # TODO: get CRS from description
    msh = {
        "coords": coords,
        "triangles": triangles,
        "quads": quads,
        "values": values,
        "description": description,
    }
    if 'crs' in grd:
        msh.update({"crs": grd["crs"]})
    if 'boundaries' in grd.keys():
        msh.update({'boundaries': grd['boundaries']})
    return msh
