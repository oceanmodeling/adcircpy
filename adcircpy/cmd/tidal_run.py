from datetime import datetime, timedelta
import logging

from pytz import timezone

from adcircpy.cmd import argument_parser
from adcircpy.cmd.basecmd import AdcircCommand


class TidalRunCommand(AdcircCommand):
    """CLI wrapper for AdcircCommand to generate tidal only runs"""

    def __init__(self, args):
        super().__init__(args)
        self.start_date = datetime.strptime(self.args.start_date, '%Y-%m-%dT%H:%M:%S')
        self.end_date = self.start_date + timedelta(days=self.args.run_days)


def main():
    args = argument_parser.get_parser('tidal').parse_args()
    if len(args.constituents) == 0:
        args.constituents = ['all']
    logging.basicConfig(
        level={'warning': logging.WARNING, 'info': logging.INFO, 'debug': logging.DEBUG,}[
            args.log_level
        ],
        format='[%(asctime)s] %(name)s %(levelname)s: %(message)s',
        force=True,
    )
    logging.Formatter.converter = lambda *args: datetime.now(tz=timezone('UTC')).timetuple()
    TidalRunCommand(args).run()
