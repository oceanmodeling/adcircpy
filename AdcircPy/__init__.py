try:
    import colored_traceback
    colored_traceback.add_hook(always=True)
except ModuleNotFoundError:
    pass

# required for plotting a large number of mesh-elements.
import matplotlib as mpl  # noqa: E402
mpl.rcParams['agg.path.chunksize'] = 10000

from AdcircPy import Model  # noqa: E402
from AdcircPy.AdcircPy import read_mesh, read_output  # noqa: E402


__author__ = "Jaime R. Calzada-Marrero"
__copyright__ = "Copyright 2019, Jaime R. Calzada-Marrero"
__credits__ = ["Saeed Moghimi",
               "Sergey Vinogradov",
               "Edward Myers",
               "Yuji Funakoshi",
               "NOAA/NOS/CSDL"]
__license__ = "GPL"
__version__ = "0.0.9"
__maintainer__ = "Jaime R. Calzada"
__email__ = "jaime.calzada@noaa.gov"
__status__ = "Development"

__all__ = ['read_mesh',
           'read_output',
           'Model']
