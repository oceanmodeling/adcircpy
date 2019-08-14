#!/usr/bin/env python
import argparse
import matplotlib.pyplot as plt
from adcircpy.mesh import AdcircMesh


class PlotMesh:

    def __init__(self):
        mesh = AdcircMesh.open(self.args.mesh_path)
        ax = mesh.make_plot(vmin=self.args.vmin, vmax=self.args.vmax)
        if self.args.show_elements:
            ax.triplot(mesh.mpl_tri, color='k', linewidth=0.07)
        plt.show()

    def parse_args(self):
        parser = argparse.ArgumentParser(
            description="Program to see a quick plot of an ADCIRC mesh.")
        parser.add_argument("mesh_path", help="ADCIRC mesh file path.")
        parser.add_argument("--show-elements", action="store_true",
                            default=False)
        parser.add_argument("--vmin", type=float)
        parser.add_argument("--vmax", type=float)
        return parser.parse_args()

    @property
    def args(self):
        if not hasattr(self, "__args"):
            self.__args = self.parse_args()
        return self.__args


def main():
    PlotMesh()


if __name__ == "__main__":
    main()
