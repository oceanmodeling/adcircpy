from importlib import util

import matplotlib as mpl
from pandas.plotting import register_matplotlib_converters

from adcircpy.driver import AdcircRun
from adcircpy.forcing import Tides
from adcircpy.mesh import AdcircMesh
from adcircpy.fort15 import Fort15

__all__ = [
    "AdcircMesh",
    "AdcircRun",
    "Tides",
    'Fort15'
]

mpl.rcParams['agg.path.chunksize'] = 10000
register_matplotlib_converters()

if util.find_spec("colored_traceback") is not None:
    import colored_traceback

    colored_traceback.add_hook(always=True)
