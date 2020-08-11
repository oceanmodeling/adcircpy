from datetime import datetime
from functools import lru_cache
import logging

from adcircpy.cmd import argument_parser
from adcircpy.cmd.basecmd import AdcircCommand

try:
    import colored_traceback

    colored_traceback.add_hook(always=True)
except ModuleNotFoundError:
    pass


class TidalRunCommand(AdcircCommand):
    """ CLI wrapper for AdcircCommand to generate tidal only runs """

    @property
    @lru_cache(maxsize=None)
    def _start_date(self):
        return datetime.strptime(self.args.start_date, "%Y-%m-%dT%H:%M")

    @property
    @lru_cache(maxsize=None)
    def _end_date(self):
        return datetime.strptime(self.args.end_date, "%Y-%m-%dT%H:%M")


def main():
    args = argument_parser.get_parser('tidal').parse_args()
    if len(args.constituents) == 0:
        args.constituents = ['all']
    logging.basicConfig(level=args.log_level)
    logging.getLogger("rasterio").setLevel(logging.WARNING)
    logging.getLogger("fiona").setLevel(logging.WARNING)
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("paramiko").setLevel(logging.WARNING)
    drv = TidalRunCommand(args)
    retv = drv.run()
    exit(retv)
