# import argparse
import logging

from adcircpy.cmd import argument_parser
from adcircpy.cmd.basecmd import AdcircCommand
from adcircpy.forcing.winds.best_track import BestTrackForcing

_logger = logging.getLogger(__name__)


class BestTrackRunCommand(AdcircCommand):

    def __init__(self, args):
        super().__init__(args)
        bt = BestTrackForcing(self._args.storm_id)
        if args.start_date is None and args.end_date is None:
            bt.clip_to_bbox(self.mesh.get_bbox(), self.mesh.crs)
            self._start_date = bt.start_date
            self._end_date = bt.end_date
        else:
            # TODO: make sure the setters can parse the arguments.
            self._start_date = args.start_date
            self._end_date = args.end_date

        self.mesh.add_forcing(bt)


def set_logger(log_level):
    logging.basicConfig(level=log_level)
    logging.getLogger("rasterio").setLevel(logging.WARNING)
    logging.getLogger("fiona").setLevel(logging.WARNING)
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("paramiko").setLevel(logging.WARNING)


def main():
    args = argument_parser.get_parser('best_track').parse_args()
    set_logger(args.log_level)
    # if len(args.constituents) == 0:
    #     args.constituents = ['all']
    exit(BestTrackRunCommand(args).run())
