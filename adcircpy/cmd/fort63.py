import argparse

import matplotlib.pyplot as plt

from adcircpy.cmd import diagnose
from adcircpy.outputs.fort63 import Fort63


def plot(fort63, args):
    fort63.index = args.index
    ax = fort63.tricontourf(
        vmin=args.vmin,
        vmax=args.vmax,
        cbar=True,
    )
    if args.plot_elements:
        ax.triplot(fort63.triangulation, color='k', linewidth=0.1)
    if args.diagnose is not None:
        elmax, speedmax, index = diagnose.parse(args.diagnose)
        ax.scatter(
            fort63.x[index],
            fort63.y[index],
            # c=elmax,
            marker='o',
            edgecolor='r',
            facecolor='none'
        )
        plt.show()


def animation(fort63, args):
    fort63.animation(save=args.save_path, show=args.show)


def export(fort63, args):
    fort63.index = args.index
    fort63.export(args.output, overwrite=args.overwrite)


def main():
    args = parse_args()
    fort63 = Fort63(
        args.fort63,
        # fort14=args.fort14
    )
    {
        'plot': plot,
        'animation': animation,
        'export': export
    }[args.mode](fort63, args)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Program to see a quick plot of an ADCIRC fort63 file."
    )
    parser.add_argument('fort63', help="Path to fort.63 file.")

    subparsers = parser.add_subparsers(dest='mode')
    subparsers.required = True

    # data plotting subparsers
    plot = subparsers.add_parser('plot')
    plot.add_argument('index', type=int, default=-1)
    plot.add_argument('--no-show', action='store_true')
    plot.add_argument('--save')
    plot.add_argument('--fps', type=int, default=2)
    _help = "Path to fort.14 file (required if fort63 files is not netcdf)."
    plot.add_argument('--fort14', help=_help)
    plot.add_argument('--title', help="Plot title override.")
    plot.add_argument('--vmin', type=float)
    plot.add_argument('--vmax', type=float)
    plot.add_argument('--start-index', type=int)
    plot.add_argument('--end-index', type=int)
    plot.add_argument("--plot-elements", action="store_true")
    plot.add_argument("--diagnose")
    plot.add_argument('--save-path')

    # export
    export = subparsers.add_parser('export')
    export.add_argument('output')
    export.add_argument('index', type=int)
    export.add_argument('--overwrite', action='store_true')

    return parser.parse_args()


if __name__ == '__main__':
    main()
