# import argparse
import logging
from adcircpy.model.winds import BestTrackForcing
from adcircpy.cmd.basecmd import AdcircCommand
from adcircpy.cmd import argument_parser


class BestTrackRunCommand(AdcircCommand):
    """
    """

    def __init__(self, args):
        super().__init__(args)
        self._init_best_track()

    def _init_best_track(self):
        bt = BestTrackForcing(self.args.storm_id)
        if self.args.start_date is None:
            bt.clip_to_bbox(self.mesh.get_bbox("EPSG:3395"))
            self.__start_date = bt.start_date
            self.__end_date = bt.end_date
        else:
            raise NotImplementedError("add custom date times?")
        self.__wind_forcing = bt

    @property
    def _wind_forcing(self):
        return self.__wind_forcing

    @property
    def _start_date(self):
        return self.__start_date

    @property
    def _end_date(self):
        return self.__end_date


def main():

    args = argument_parser.get_parser('best_track').parse_args()
    if len(args.constituents) == 0:
        args.constituents = ['all']
    logging.basicConfig(level=args.log_level)
    logging.getLogger("rasterio").setLevel(logging.WARNING)
    logging.getLogger("fiona").setLevel(logging.WARNING)
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("paramiko").setLevel(logging.WARNING)
    drv = BestTrackRunCommand(args)
    retv = drv.run()
    exit(retv)

    # args = parse_args()
    # drv = BestTrackRunCommand(args)
    # retv = drv.run()
    # exit(retv)
