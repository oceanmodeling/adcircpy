# import argparse
import logging

from adcircpy.cmd import argument_parser
from adcircpy.cmd.basecmd import AdcircCommand
from adcircpy.forcing.winds.best_track import BestTrackForcing


class BestTrackRunCommand(AdcircCommand):

    def __init__(self, args):
        super().__init__(args)
        bt = BestTrackForcing(self._args.storm_id)
        self.mesh.add_forcing(bt)
        if self._args.start_date is None:
            self._start_date = bt.start_date
            self._end_date = bt.end_date
        else:
            raise NotImplementedError("add custom date times?")


def main():
    args = argument_parser.get_parser('best_track').parse_args()
    # if len(args.constituents) == 0:
    #     args.constituents = ['all']
    logging.basicConfig(level=args.log_level)
    logging.getLogger("rasterio").setLevel(logging.WARNING)
    logging.getLogger("fiona").setLevel(logging.WARNING)
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("paramiko").setLevel(logging.WARNING)
    exit(BestTrackRunCommand(args).run())
