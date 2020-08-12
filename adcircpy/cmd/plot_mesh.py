#!/usr/bin/env python
import argparse
import pathlib

import matplotlib.pyplot as plt

from adcircpy import AdcircMesh


class PlotMeshCommand:

    def __init__(self, args):
        self._args = args

    def run(self):
        self._make_plot()
        self._make_triplot()
        self._make_boundary_plot()
        plt.gca().axis('scaled')
        plt.show()
        return 0

    def _make_plot(self):
        if not self.args.no_topobathy:
            self.mesh.make_plot(
                axes=self.ax,
                vmin=self.args.vmin,
                vmax=self.args.vmax,
                # levels=self.args.levels
            )
        if self.args.diagnose is not None:
            elmax, speedmax, index = diagnose(self.args.diagnose)
            self.ax.scatter(
                self.mesh.x[index],
                self.mesh.y[index],
                # c=elmax,
                marker='o',
                edgecolor='r',
                facecolor='none'
            )

    def _make_triplot(self):
        if self.args.show_elements:
            self.ax.triplot(self.mesh.triangulation, color='k', linewidth=0.07)

    def _make_boundary_plot(self):
        if self.args.plot_boundaries:
            self.mesh.plot_ocean_boundaries(axes=self.ax)

    @property
    def args(self):
        return self._args

    @property
    def mesh(self):
        try:
            return self.__mesh
        except AttributeError:
            self.__mesh = AdcircMesh.open(self.args.mesh)
            return self.__mesh

    @property
    def fig(self):
        try:
            return self.__fig
        except AttributeError:
            self.__fig = plt.figure()
            return self.__fig

    @property
    def ax(self):
        try:
            return self.__ax
        except AttributeError:
            self.__ax = self.fig.add_subplot(111)
            return self.__ax

    @property
    def _args(self):
        return self.__args

    @_args.setter
    def _args(self, args):
        self.__args = args


def diagnose(logfile):
    import numpy as np
    logfile = pathlib.Path(logfile).resolve()
    with open(logfile, 'r') as f:
        lines = "".join(f.readlines())
    elmax = list()
    speedmax = list()
    index = list()
    _lines = lines.split(
        '** ERROR: Elevation.gt.ErrorElev, ADCIRC stopping. **\n')
    line0 = "".join(_lines[0]).split('\n')
    for line in line0:
        if "** WARNING: Elevation.gt.WarnElev **" in line:
            elmax.append(float(line.split("AT NODE")[0].split("=")[-1]))
            speedmax.append(
                float(line.split("SPEEDMAX =")[0].split("AT NODE")[-1]))
            index.append(np.abs(int(
                line.split("AT NODE")[-1].split("ON MYPROC")[0].strip())) - 1)
    return elmax, speedmax, index


def parse_args():
    parser = argparse.ArgumentParser(
        description="Program to see a quick plot of an ADCIRC mesh.")
    parser.add_argument("mesh", help="ADCIRC mesh file path.")
    parser.add_argument("--show-elements", action="store_true",
                        default=False)
    parser.add_argument("--no-topobathy", action="store_true",
                        default=False)
    parser.add_argument("--vmin", type=float)
    parser.add_argument("--vmax", type=float)
    parser.add_argument("--plot-boundaries", action="store_true")
    parser.add_argument("--diagnose")
    return parser.parse_args()


def main():
    exit(PlotMeshCommand(parse_args()).run())


if __name__ == "__main__":
    main()
