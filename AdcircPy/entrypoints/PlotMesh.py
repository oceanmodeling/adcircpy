import argparse
from AdcircPy import AdcircPy


class PlotMesh(object):
    def __init__(self):
        self.parse_args()
        self.read_mesh()
        self.plot_mesh()

    def parse_args(self):
        parser = argparse.ArgumentParser(
                        description="Program to see a quick plot of an "
                                    + "ADCIRC mesh.")
        parser.add_argument("mesh_path", help="ADCIRC mesh file path.")
        parser.add_argument("--show-elements", action="store_true",
                            default=False)
        parser.add_argument("--vmin", type=float)
        parser.add_argument("--vmax", type=float)
        self.args = parser.parse_args()

    def read_mesh(self):
        self.fort14 = AdcircPy.read_mesh(self.args.mesh_path)

    def plot_mesh(self):
        self.fort14.make_plot(show=True, elements=self.args.show_elements,
                              vmin=self.args.vmin, vmax=self.args.vmax)


def main():
    PlotMesh()


if __name__ == "__main__":
    main()
