import sys

import gartersnake
from setuptools import config, find_packages, setup

DEPENDENCIES = {
    'appdirs': [],
    'geopandas': ['gdal', 'fiona'],
    'haversine': [],
    'matplotlib': [],
    'netCDF4': [],
    'numpy': [],
    'pandas': [],
    'paramiko': [],
    'psutil': [],
    'pyproj>=2.6': [],
    'requests': [],
    'scipy': [],
    'shapely': [],
    'stormevents>=1.2': [],
    'utm': [],
}

MISSING_DEPENDENCIES = gartersnake.missing_requirements(DEPENDENCIES)

if len(MISSING_DEPENDENCIES) > 0:
    print(f'{len(MISSING_DEPENDENCIES)} (out of {len(DEPENDENCIES)}) dependencies are missing')

if len(MISSING_DEPENDENCIES) > 0 and gartersnake.is_conda():
    print(f'found conda environment at {sys.prefix}')
    gartersnake.install_conda_requirements(MISSING_DEPENDENCIES)
    MISSING_DEPENDENCIES = gartersnake.missing_requirements(DEPENDENCIES)

if len(MISSING_DEPENDENCIES) > 0 and gartersnake.is_windows():
    gartersnake.install_windows_requirements(MISSING_DEPENDENCIES)
    MISSING_DEPENDENCIES = gartersnake.missing_requirements(DEPENDENCIES)

__version__ = gartersnake.vcs_version()
print(f'using version {__version__}')

metadata = config.read_configuration('setup.cfg')['metadata']

setup(
    **metadata,
    version=__version__,
    packages=find_packages(),
    python_requires='>=3.6',
    setup_requires=['dunamai', 'gartersnake', 'setuptools>=41.2'],
    install_requires=list(DEPENDENCIES),
    # test and development dependencies
    extras_require={
        'testing': [
            'pooch',
            'pytest',
            'pytest-cov',
            'pytest-mock',
            'pytest-socket',
            'pytest-xdist',
        ],
        'development': ['dunamai', 'flake8', 'isort', 'oitnb'],
        'documentation': [
            'dunamai',
            'm2r2',
            'sphinx',
            'sphinx-rtd-theme',
            'sphinxcontrib-programoutput',
            'sphinxcontrib-bibtex',
        ],
    },
    entry_points={
        'console_scripts': [
            'tidal_run=adcircpy.cmd.tidal_run:main',
            'best_track_run=adcircpy.cmd.best_track_run:main',
            'best_track_file=adcircpy.cmd.best_track_file:main',
            'plot_mesh=adcircpy.cmd.plot_mesh:main',
            'plot_maxele=adcircpy.cmd.plot_maxele:main',
            'plot_fort61=adcircpy.cmd.plot_fort61:main',
            'fort63=adcircpy.cmd.fort63:main',
            'tide_gen=adcircpy.cmd.tide_gen:main',
        ]
    },
)
