from os import PathLike
from pathlib import Path

from adcircpy.forcing.waves import WaveForcing

import warnings

from adcircpy.warnings import warn_adcirc, UnsupportedFeatureWarning

class SWANForcing(WaveForcing):
    def __init__(self, nrs: int = 3, \
            interval_seconds: int = 600):
        super().__init__(nrs=nrs, interval_seconds=interval_seconds)

    def write(self, directory: PathLike, overwrite: bool = False):
        txt=f'ADCIRCPy does not currently support writing '\
                +f'coupled ADCIRC+SWAN input files (fort.26 and swaninit). '\
                +f'User will have to create these files themselves to run '\
                +f'padcswan.'
        warn_adcirc(txt, UnsupportedFeatureWarning)
        pass
