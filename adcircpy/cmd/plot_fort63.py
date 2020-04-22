import argparse
import pathlib
import matplotlib.pyplot as plt
from adcircpy.cmd import diagnose
from adcircpy.outputs.fort63 import Fort63


def main():
    args = parse_args()
    fort63 = Fort63.open(args.fort63, fort14=args.fort14)

    if args.index is not None:
        ax = fort63.make_plot(
            args.index,
            title=args.title,
            vmin=args.vmin,
            vmax=args.vmax,
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

    else:
        anim = fort63.make_movie()
        if args.save_path is not None:
            anim.save(
                pathlib.Path(args.save_path).resolve(),
                writer='imagemagick',
                # fps=60
                )
        else:
            plt.show()


def parse_args():
    parser = argparse.ArgumentParser(
        description="Program to see a quick plot of an ADCIRC fort63 file."
        )
    parser.add_argument(
        'fort63',
        help="Path to fort.63 file."
        )
    parser.add_argument(
        '--index', type=int,
        help='')
    parser.add_argument(
        '--fort14',
        help="Path to fort.14 file (required if fort63 files is not netcdf)."
        )
    parser.add_argument(
        '--title',
        help="Plot title override."
        )
    parser.add_argument(
        '--vmin',
        type=float
        )
    parser.add_argument(
        '--vmax',
        type=float
        )

    parser.add_argument(
        '--start-index',
        type=float
        )
    parser.add_argument(
        '--end-index',
        type=float
        )
    parser.add_argument("--plot-elements", action="store_true")
    parser.add_argument("--diagnose")
    parser.add_argument('--save-path')
    return parser.parse_args()


if __name__ == '__main__':
    main()
