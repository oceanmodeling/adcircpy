# import argparse
from datetime import datetime, timedelta
import logging

from pytz import timezone

from adcircpy.cmd import argument_parser
from adcircpy.cmd.basecmd import AdcircCommand
from adcircpy.forcing.winds.best_track import BestTrackForcing
from adcircpy.utilities import get_logger

LOGGER = get_logger(__name__)


class BestTrackRunCommand(AdcircCommand):
    def __init__(self, args):

        LOGGER.info('Init BestTrackRunCommand')
        super().__init__(args)

        LOGGER.info(f'Init BestTrackForcing for {self.args.storm_id}')
        bt = BestTrackForcing(self.args.storm_id)

        LOGGER.info('Clip BestTrackForcing to bbox')
        if self.args.clip:
            bt.clip_to_bbox(self.mesh.get_bbox(output_type='bbox'), self.mesh.crs)

        if args.start_date is None:
            self.start_date = bt.start_date
        else:
            self.start_date = datetime.strptime(args.start_date, '%%Y-%%m-%%dT%%H')

        if args.run_days is None:
            self.end_date = bt.end_date
        else:
            self.end_date = self.start_date + timedelta(days=args.run_days)

        bt.start_date = self.start_date
        bt.end_date = self.end_date

        self.mesh.add_forcing(bt)


def main():
    args = argument_parser.get_parser('best_track').parse_args()
    logging.basicConfig(
        level={'warning': logging.WARNING, 'info': logging.INFO, 'debug': logging.DEBUG,}[
            args.log_level
        ],
        format='[%(asctime)s] %(name)s %(levelname)s: %(message)s',
        # force=True,
    )
    logging.Formatter.converter = lambda *args: datetime.now(tz=timezone('UTC')).timetuple()
    BestTrackRunCommand(args).run()
