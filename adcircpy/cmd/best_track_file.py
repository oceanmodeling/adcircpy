#! /usr/bin/env python
import argparse

from adcircpy.forcing.winds.best_track import BestTrackForcing


def parse_args():
    parser = argparse.ArgumentParser(
        description="Program to generate fort.22 from HURDAT2.")
    parser.add_argument(
        "storm_id",
        help="Storm id from HURDAT2 database.")
    parser.add_argument(
        "--save-path",
        help="Output path for fort.22 file.")
    parser.add_argument(
        '--start-date', help="format is %Y%m%d%H")
    parser.add_argument(
        '--end-date', help="format is %Y%m%d%H")
    parser.add_argument(
        '--quiet', '-q', action='store_true', default=False)
    parser.add_argument(
        '--six-hourly', '-s', action='store_true', default=False,
        help='Filter data to six-hourly only.')
    parser.add_argument(
        '--show-HU-only', '-hu', action='store_true', default=False,
        help='Removes both tropical storm (TS) and extratropical (EX) entries.'
    )
    parser.add_argument(
        '--remove-TS', action='store_true', default=False,
        help='Removes tropical storm (TS) entries.')
    parser.add_argument(
        '--remove-EX', action='store_true', default=False,
        help='Removes extratropical (EX) entries.')
    parser.add_argument(
        '--plot-track', action='store_true', default=False,
        help='Shows a simple plot of the track.')
    return parser.parse_args()


def main():
    args = parse_args()
    bt = BestTrackForcing(args.storm_id)
    print(bt.fort22)
    # # set custom start date
    # if args.start_date is not None:
    #     bt.start_date = args.start_date
    # # set custom end date
    # if args.end_date is not None:
    #     bt.end_date = args.end_date
    # # filter non-six-hourly data
    # if args.six_hourly:
    #     bt.remove_non_six_hourly()
    # # remove non-HU data
    # if args.show_HU_only:
    #     bt.only_HU()
    # # remove TS data
    # if args.remove_TS:
    #     bt.remove_TS()
    # # remove EX data
    # if args.remove_EX:
    #     bt.remove_EX()
    # # print fort22
    # if args.quiet is False:
    #     print(bt.fort22)
    # # show cheap plot
    # if args.plot_track:
    #     bt.plot_track()
    # # save fort22
    # if args.save_path is not None:
    #     bt.dump(args.save_path)
    #     if args.quiet is False:
    #         print('File written to: {}'.format(args.save_path))


if __name__ == "__main__":
    main()
