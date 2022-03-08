from datetime import datetime
from enum import Enum
import logging
from os import PathLike
import pathlib


class MeshGeometryType(Enum):
    TRIANGLE = 'E3T'
    QUADRILATERAL = 'E4Q'
    HEXAGON = 'E6T'
    OCTAGON = 'E8Q'
    NONAGON = 'E9Q'


def read(path):
    mesh = {}
    with open(pathlib.Path(path), 'r') as f:
        f.readline()

        for line in f.readlines():
            line = line.split()

            geom_type = line[0]
            if geom_type not in mesh:
                mesh[geom_type] = {}

            if geom_type in ['E3T', 'E4Q']:
                mesh[geom_type].update({line[1]: line[2:]})
            elif geom_type in ['ND']:
                mesh[geom_type].update(
                    {line[1]: (list(map(float, line[2:-1])), float(line[-1]))}
                )
    return mesh


def write(mesh: {str: {str: (float, float)}}, path: PathLike, overwrite: bool = False):
    if not isinstance(path, pathlib.Path):
        path = pathlib.Path(path)

    triangles = mesh[MeshGeometryType.TRIANGLE.value]
    triangles.insert(0, 'type', MeshGeometryType.TRIANGLE.value)
    triangles.insert(1, 'id', triangles.index)
    quadrilaterals = mesh[MeshGeometryType.QUADRILATERAL.value]
    quadrilaterals.insert(0, 'type', MeshGeometryType.QUADRILATERAL.value)
    quadrilaterals.insert(1, 'id', quadrilaterals.index)
    nodes = mesh['ND']
    nodes.insert(0, 'type', 'ND')
    nodes.insert(1, 'id', nodes.index)

    if 'boundaries' in mesh:
        boundaries = mesh['boundaries']
        boundaries.insert(0, 'type', 'NS')
        boundaries.iloc[:, 2:] *= -1
    else:
        boundaries = None

    def float_format(value: float):
        return f'{value:<.16E}'

    if overwrite or not path.exists():
        with open(path, 'w') as f:
            f.write('MESH2D\n')

            if len(triangles) > 0:
                logging.debug('writing triangles')
                start_time = datetime.now()
                triangles.to_string(f, header=False, index=False, justify='left')
                f.write('\n')
                logging.debug(f'wrote triangles in {datetime.now() - start_time}')

            if len(quadrilaterals) > 0:
                logging.debug('writing quadrilaterals')
                start_time = datetime.now()
                quadrilaterals.to_string(f, header=False, index=False, justify='left')
                f.write('\n')
                logging.debug(f'wrote quadrilaterals in {datetime.now() - start_time}')

            logging.debug('writing nodes')
            start_time = datetime.now()
            nodes.to_string(
                f, header=False, index=False, justify='left', float_format=float_format
            )
            f.write('\n')
            logging.debug(f'wrote nodes in {datetime.now() - start_time}')

            if boundaries in mesh:
                logging.debug('writing boundaries')
                start_time = datetime.now()
                boundaries.to_string(f, header=False, index=False, justify='left')
                f.write('\n')
                logging.debug(f'wrote boundaries in {datetime.now() - start_time}')

        return 0  # for unittests
    else:
        logging.debug(f'skipping existing file "{path}"')
        return 1
