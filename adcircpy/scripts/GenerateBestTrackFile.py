#! /usr/bin/env python
import argparse
from AdcircPy.model import BestTrackForcing


class GenerateBestTrackFile(object):

    def __init__(self):
        # fetch storm data from internet
        bt = BestTrackForcing(self.args.storm_id)
        # set custom start date
        if self.args.start_date is not None:
            bt.start_date = self.args.start_date
        # set custom end date
        if self.args.end_date is not None:
            bt.end_date = self.args.end_date
        # filter non-six-hourly data
        if self.args.six_hourly:
            bt.remove_non_six_hourly()
        # remove non-HU data
        if self.args.show_HU_only:
            bt.only_HU()
        # remove TS data
        if self.args.remove_TS:
            bt.remove_TS()
        # remove EX data
        if self.args.remove_EX:
            bt.remove_EX()
        # print fort22
        if self.args.quiet is False:
            print(bt.fort22)
        # show cheap plot
        if self.args.plot_track:
            bt.plot_track()
        # save fort22
        if self.args.save_path is not None:
            bt.dump(self.args.save_path)
            if self.args.quiet is False:
                print('File written to: {}'.format(self.args.save_path))

    @staticmethod
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
            help='Filter data to HU entries only.')
        parser.add_argument(
            '--remove-TS', action='store_true', default=False,
            help='Removes TS data from dataset.')
        parser.add_argument(
            '--remove-EX', action='store_true', default=False,
            help='Shows a cheap plot of the track.')
        parser.add_argument(
            '--plot-track', action='store_true', default=False,
            help='Shows a cheap plot of the track.')
        return parser.parse_args()

    @property
    def args(self):
        if not hasattr(self, "__args"):
            self.__args = self.parse_args()
        return self.__args


def main():
    GenerateBestTrackFile()


if __name__ == "__main__":
    main()
