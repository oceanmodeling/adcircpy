#!/usr/bin/env python
import argparse
import pathlib

import matplotlib.pyplot as plt

from adcircpy import AdcircMesh


class PlotMeshCommand:
    def __init__(self, args: argparse.Namespace):
        mesh = AdcircMesh.open(args.mesh, crs=args.crs)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if args.no_topobathy is False:
            mesh.make_plot(axes=ax, vmin=args.vmin, vmax=args.vmax)
        if args.show_elements:
            mesh.triplot(axes=ax)
        if args.plot_boundaries:
            mesh.boundaries.gdf.plot(ax=ax)
        plt.show(block=True)


def diagnose(logfile):
    import numpy as np

    logfile = pathlib.Path(logfile).resolve()
    with open(logfile, 'r') as f:
        lines = "".join(f.readlines())
    elmax = list()
    speedmax = list()
    index = list()
    _lines = lines.split('** ERROR: Elevation.gt.ErrorElev, ADCIRC stopping. **\n')
    line0 = "".join(_lines[0]).split('\n')
    for line in line0:
        if '** WARNING: Elevation.gt.WarnElev **' in line:
            elmax.append(float(line.split('AT NODE')[0].split('=')[-1]))
            speedmax.append(float(line.split('SPEEDMAX =')[0].split('AT NODE')[-1]))
            index.append(
                np.abs(int(line.split('AT NODE')[-1].split('ON MYPROC')[0].strip())) - 1
            )
    return elmax, speedmax, index


def parse_args():
    parser = argparse.ArgumentParser(
        description='Program to see a quick plot of an ADCIRC mesh.'
    )
    parser.add_argument('mesh', help='ADCIRC mesh file path.')
    parser.add_argument('--crs', help='ADCIRC mesh crs.')
    parser.add_argument('--show-elements', action='store_true', default=False)
    parser.add_argument('--no-topobathy', action='store_true', default=False)
    parser.add_argument('--vmin', type=float)
    parser.add_argument('--vmax', type=float)
    parser.add_argument('--plot-boundaries', action='store_true')
    parser.add_argument('--diagnose')
    return parser.parse_args()


def main():
    PlotMeshCommand(parse_args())


if __name__ == '__main__':
    main()
