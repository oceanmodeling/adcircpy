#! /usr/bin/env python
import argparse
from AdcircPy.Winds import BestTrackForcing
# from AdcircPy.Winds import HURDAT2


class GenerateBestTrackFile(object):
    def __init__(self):
        self._parse_args()
        self._get_fort22()
        self._print_fort22()
        self._save_fort22()

    def _parse_args(self):
        parser = argparse.ArgumentParser(description="Program to generate "
                                         + "fort.22 from HURDAT2.")
        parser.add_argument("storm_id",  help="Storm id from HURDAT2 database."
                            )
        parser.add_argument("--save-path", help="Output path for fort.22 file."
                            )
        parser.add_argument('--verbose', '-v', dest='verbose',
                            action='store_true')
        parser.add_argument('--quiet', dest='verbose', action='store_false')
        parser.set_defaults(verbose=False)
        self.args = parser.parse_args()

    def _get_fort22(self):
        # self.fort22 = HURDAT2(self.args.storm_id)
        # print(self.fort22.)
        self.fort22 = BestTrackForcing(self.args.storm_id)

    def _print_fort22(self):
        if self.args.verbose is True:
            print(self.fort22.get_fort22())

    def _save_fort22(self):
        if self.args.save_path is not None:
            self.fort22.dump(self.args.save_path)
            if self.args.verbose is True:
                print('File written to: {}'.format(self.args.save_path))


def main():
    GenerateBestTrackFile()


if __name__ == "__main__":
    main()
