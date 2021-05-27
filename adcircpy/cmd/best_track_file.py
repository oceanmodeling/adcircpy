#! /usr/bin/env python
import argparse
from datetime import datetime
from pathlib import Path

from adcircpy.forcing.winds.best_track import BestTrackForcing


def parse_args():
    parser = argparse.ArgumentParser(
        description='generate `fort.22` information from HURDAT2 data'
    )
    parser.add_argument(
        'storm',
        help='Can be StormYYYY (eg. Sandy2012), '
        'a storm id from the HURDAT2 table (ftp://ftp.nhc.noaa.gov/atcf/archive/storm.table), or'
        'a file path to a track file in HURDAT2 format',
    )
    parser.add_argument('--save-path', help='path to which to write fort.22')
    parser.add_argument('--start-date', help='format is %%Y%%m%%d%%H')
    parser.add_argument('--end-date', help='format is %%Y%%m%%d%%H')
    parser.add_argument(
        '--quiet', '-q', action='store_true', default=False, help='suppress console output',
    )
    parser.add_argument(
        '--plot-track',
        action='store_true',
        default=False,
        help='show a simple plot of the track',
    )
    parser.add_argument('--nws', help='which NWS forcing option to use')
    return parser.parse_args()


def main():
    args = parse_args()

    if 'fort.22' in args.storm:
        bt = BestTrackForcing.from_fort22(
            args.storm, nws=int(args.nws), start_date=args.start_date, end_date=args.end_date,
        )
    else:
        bt = BestTrackForcing(
            args.storm,
            nws=int(args.nws) if args.nws is not None else 20,
            start_date=datetime.strptime(args.start_date, '%Y%m%d%H')
            if args.start_date is not None
            else None,
            end_date=datetime.strptime(args.end_date, '%Y%m%d%H')
            if args.end_date is not None
            else None,
        )

    # print fort22
    if not args.quiet:
        print(str(bt))

    # show cheap plot
    if args.plot_track:
        bt.plot_track(show=True)

    # save fort22
    if args.save_path is not None:
        with open(Path(args.save_path), 'w') as output_file:
            output_file.write(str(bt))


if __name__ == '__main__':
    main()
