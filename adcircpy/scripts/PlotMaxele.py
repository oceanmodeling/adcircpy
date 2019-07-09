import argparse
import matplotlib.pyplot as plt
from AdcircPy import read_output


class PlotMaxele(object):

    def __init__(self, args):
        maxele = read_output(args.maxele)
        maxele.make_plot(title=args.title,
                         vmin=args.vmin,
                         vmax=args.vmax)
        plt.show()


def get_parser():
    parser = argparse.ArgumentParser(
        description="Program to see a quick plot of an ADCIRC maxele file.")
    parser.add_argument(
        'maxele',
        help="Path to maxele file.")
    parser.add_argument(
        '--fort14',
        help="Path to fort.14 file (required if maxele files is not netcdf).")
    parser.add_argument('--title', help="Plot title override.")
    parser.add_argument('--vmin', type=float)
    parser.add_argument('--vmax', type=float)
    return parser


def main():
    PlotMaxele(get_parser().parse_args())
