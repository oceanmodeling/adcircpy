import matplotlib as mpl
from pandas.plotting import register_matplotlib_converters

from adcircpy.driver import AdcircRun
from adcircpy.forcing import TidalSource, Tides, WaveForcing, WindForcing
from adcircpy.fort15 import Fort15
from adcircpy.mesh import AdcircMesh

__all__ = [
    'AdcircMesh',
    'AdcircRun',
    'Tides',
    'TidalSource',
    'WaveForcing',
    'WindForcing',
    'Fort15',
]

mpl.rcParams['agg.path.chunksize'] = 10000
register_matplotlib_converters()

try:
    import colored_traceback

    colored_traceback.add_hook(always=True)
except ImportError:
    pass
