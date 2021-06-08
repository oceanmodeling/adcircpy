import argparse

import matplotlib.pyplot as plt

from adcircpy.outputs import Maxele, MaximumElevationTimes


def parse_args():
    parser = argparse.ArgumentParser(
        description='Program to generate a plot of an ADCIRC `maxele` file.'
    )
    parser.add_argument('maxele', help='Path to maxele file.')
    parser.add_argument(
        '--plot_timestep_of_maxele',
        help='Plot the timestep of the maximum elevation value',
        action='store_true',
    )
    parser.add_argument(
        '--fort14', help='Path to fort.14 file (required if maxele files is not netcdf).',
    )
    parser.add_argument('--title', help='Plot title override.')
    parser.add_argument('--vmin', type=float)
    parser.add_argument('--vmax', type=float)
    parser.add_argument('--cmap', type=str, default='jet')
    parser.add_argument('--levels', type=int, default=256)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.plot_timestep_of_maxele:
        maxele = MaximumElevationTimes(args.maxele)
    else:
        maxele = Maxele(args.maxele)
    maxele.tricontourf(
        vmin=args.vmin, vmax=args.vmax, cmap=args.cmap, levels=args.levels, cbar=True
    )
    plt.show()
