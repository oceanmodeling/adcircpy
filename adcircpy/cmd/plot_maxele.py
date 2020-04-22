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
    try:
        maxele = Maxele.from_netcdf(args.maxele)
    except OSError:
        if args.fort14 is None:
            raise IOError('Must pass --fort14 when plotting ascii files.')
        maxele = Maxele.from_ascii(args.maxele, args.fort14)
    maxele.make_plot(title=args.title,
                     vmin=args.vmin,
                     vmax=args.vmax)
    plt.show()
