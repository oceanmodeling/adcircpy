#! /usr/bin/env python
import argparse
import os
from pathlib import Path

import matplotlib.pyplot as plt
import pandas
from stormevents.coops import COOPS_Station

from adcircpy.outputs import Fort61


def parse_args():
    parser = argparse.ArgumentParser(
        description='Program to see a quick plot of an ADCIRC mesh.'
    )
    parser.add_argument('path', help='Path to ADCIRC fort.61 or fort.61.nc file.')
    parser.add_argument(
        'vertical_datum',
        choices=['MHHW', 'MHW', 'MTL', 'MSL', 'MLW', 'MLLW', 'NAVD88', 'STND'],
        help='Tidal station datum, must match vertical datum of mesh.',
    )
    show = parser.add_mutually_exclusive_group(required=False)
    show.add_argument(
        '--show',
        dest='show',
        action='store_true',
        help='Shows plots to screen as they are generated (default).',
    )
    show.add_argument(
        '--no-show',
        dest='show',
        action='store_false',
        help='Prevents the plots from showing to screen. '
        + 'Useful for only saving the plots without showing them.',
    )
    parser.add_argument('--coops-only', action='store_true', help='coops plots to screen.')
    parser.add_argument(
        '--save',
        help='Directory where to save plots. '
        + "Will be created if it doesn't exist. "
        + 'It will also overwrite files unles --resume-save is used.',
    )
    parser.add_argument(
        '--resume-save',
        action='store_true',
        help='Directory where to save plots. ' + "Will be created if it doesn't exist.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    fort61 = Fort61(args.path)

    start_date = fort61.datetime[0]
    end_date = fort61.datetime[-1]

    for station_id, data in fort61:
        if args.save is not None:
            fname = str(Path(str(Path(args.save)) + '/{}.png'.format(station_id)))
            if args.resume_save and os.path.isfile(fname):
                continue

        station = COOPS_Station(station_id)

        station_data = station.get(
            start_date=start_date, end_date=end_date, datum=args.vertical_datum
        )

        if args.coops_only and not pandas.isna(station_data['v']).all():
            plt.plot(
                station_data['t'], station_data['v'], label='COOPS', color='b', linewidth=0.7
            )
            plt.plot(
                fort61.datetime, data['values'], label='ADCIRC', color='r', linewidth=0.7,
            )
        elif not args.coops_only:
            plt.plot(
                fort61.datetime, data['values'], label='ADCIRC', color='r', linewidth=0.7,
            )

        fig = plt.gcf()
        if fig.get_axes():
            fig.set_size_inches(18.5, 10.5)
            fig.gca().set_xlim(fort61.datetime[0], fort61.datetime[-1])
            plt.ylabel('water level [meters, {}]'.format(args.vertical_datum))
            plt.title('{}\n{}'.format(station.nos_id, station.name))
            plt.legend()
            if args.save is not None:
                os.makedirs(str(Path(args.save)), exist_ok=True)
                fig.savefig(fname, dpi=300, bbox_inches='tight')
            if args.show:
                plt.show()
        plt.close(fig)


if __name__ == '__main__':
    main()
