from importlib import util
from pandas.plotting import register_matplotlib_converters
import matplotlib as mpl
from adcircpy.mesh import AdcircMesh
from adcircpy.driver import AdcircRun
from adcircpy.forcing import Tides


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
