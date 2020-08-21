from importlib import util
import logging
import os
import sys

import matplotlib as mpl
from pandas.plotting import register_matplotlib_converters

from adcircpy.driver import AdcircRun
from adcircpy.forcing import Tides
from adcircpy.mesh import AdcircMesh

__all__ = [
    "AdcircMesh",
    "AdcircRun",
    "Tides",
]

mpl.rcParams['agg.path.chunksize'] = 10000
register_matplotlib_converters()

if util.find_spec("colored_traceback") is not None:
    import colored_traceback

    colored_traceback.add_hook(always=True)

from ._version import get_versions  # noqa: E402

__version__ = get_versions()['version']
del get_versions


def repository_root(path: str = None) -> str:
    if path is None:
        path = __file__
    if os.path.isfile(path):
        path = os.path.dirname(path)
    if '.git' in os.listdir(path):
        return path
    else:
        return repository_root(os.path.dirname(path))


def get_logger(
        name: str,
        log_filename: str = None,
        file_level: int = None,
        console_level: int = None,
        log_format: str = None,
        name_field_length: int = 15
) -> logging.Logger:
    if file_level is None:
        file_level = logging.DEBUG
    if console_level is None:
        console_level = logging.INFO
    logger = logging.getLogger(name)

    # check if logger is already configured
    if logger.level == logging.NOTSET and len(logger.handlers) == 0:
        # check if logger has a parent
        if '.' in name:
            logger.parent = get_logger(name.rsplit('.', 1)[0])
        else:
            # otherwise create a new split-console logger
            logger.setLevel(logging.DEBUG)
            if console_level != logging.NOTSET:
                if console_level <= logging.INFO:
                    class LoggingOutputFilter(logging.Filter):
                        def filter(self, rec):
                            return rec.levelno in (logging.DEBUG, logging.INFO)

                    console_output = logging.StreamHandler(sys.stdout)
                    console_output.setLevel(console_level)
                    console_output.addFilter(LoggingOutputFilter())
                    logger.addHandler(console_output)

                console_errors = logging.StreamHandler(sys.stderr)
                console_errors.setLevel(max((console_level, logging.WARNING)))
                logger.addHandler(console_errors)

    if log_filename is not None:
        file_handler = logging.FileHandler(log_filename)
        file_handler.setLevel(file_level)
        for existing_file_handler in [handler for handler in logger.handlers if
                                      type(handler) is logging.FileHandler]:
            logger.removeHandler(existing_file_handler)
        logger.addHandler(file_handler)

    if log_format is None:
        log_format = f'[%(asctime)s] %(name)-{name_field_length}s %(levelname)-8s: %(message)s'
    log_formatter = logging.Formatter(log_format)
    for handler in logger.handlers:
        handler.setFormatter(log_formatter)

    return logger
