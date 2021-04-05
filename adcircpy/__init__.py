from importlib import util

import matplotlib as mpl
from pandas.plotting import register_matplotlib_converters

from .driver import AdcircRun, Fort15
from .forcing import *
from .mesh import AdcircMesh

mpl.rcParams['agg.path.chunksize'] = 10000
register_matplotlib_converters()

if util.find_spec("colored_traceback") is not None:
    import colored_traceback
    colored_traceback.add_hook(always=True)
