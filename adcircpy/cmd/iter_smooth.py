#!/usr/bin/env python
"""
Takes as input the logfile of an ADCIRC run and attempts to avoid the blowup
by assigning the minimum topobathy value of the group of nodes where the blowup
ocurred.
"""
import argparse
from glob import glob
import pathlib

import matplotlib.pyplot as plt

from adcircpy.cmd import diagnose
from adcircpy.mesh import AdcircMesh


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('base_dir')
    parser.add_argument('--log-filename', default='sbatch.log')
    args = parser.parse_args()
    args.base_dir = pathlib.Path(args.base_dir).resolve()
    assert args.base_dir.is_dir()
    return args


def main():
    args = parse_args()
    iter_dir = args.base_dir / 'iter'
    iter_dir.mkdir(exist_ok=True)
    logs = glob(str(args.base_dir / '**' / args.log_filename), recursive=True)
    base_dir = args.base_dir
    if len(logs) == 0:
        msg = 'No log file found!'
        raise Exception(msg)
    elif len(logs) == 1:
        log_file = base_dir / args.log_filename
    else:
        msg = 'More than 1 logfile'
        raise NotImplementedError(msg)

    elmax, speedmax, indexes = diagnose.parse(log_file)

    if len(indexes) == 0:
        msg = 'Congratulations, your mesh did not blowup with ADCIRC. '
        msg += " That's a feat."
        print(msg)
        exit()
    mesh = AdcircMesh.open(base_dir / 'fort.14')
    # mesh.values[indexes] = np.min(mesh.values[indexes])
    ax = mesh.make_plot()
    ax.triplot(mesh.triangulation, color='k', linewidth=0.05)
    ax.scatter(mesh.x[indexes], mesh.y[indexes], edgecolor='r', facecolor='none')
    plt.show()
    return 0


if __name__ == '__main__':
    exit(main())
