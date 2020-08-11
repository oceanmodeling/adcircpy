import argparse

import matplotlib.pyplot as plt

from adcircpy.outputs import Maxele


def parse_args():
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
    return parser.parse_args()


def main():
    args = parse_args()
    maxele = Maxele(args.maxele)
    maxele.tricontourf()
    plt.show()
