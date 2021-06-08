#!/usr/bin/env python
import logging
import pathlib

import setuptools

try:
    try:
        from dunamai import Version
    except ImportError:
        import subprocess
        import sys

        subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'dunamai'])
        from dunamai import Version

    version = Version.from_any_vcs().serialize()
except RuntimeError as error:
    logging.exception(error)
    version = '0.0.0'

logging.info(f'using version {version}')

metadata = setuptools.config.read_configuration(
    pathlib.Path(__file__).parent.absolute() / 'setup.cfg'
)['metadata']
setuptools.setup(
    name=metadata['name'],
    version=version,
    author=metadata['author'],
    author_email=metadata['author_email'],
    description=metadata['description'],
    long_description=metadata['long_description'],
    long_description_content_type='text/markdown',
    url=metadata['url'],
    packages=setuptools.find_packages(),
    python_requires='>=3.6',
    setup_requires=['dunamai', 'requests', 'setuptools>=41.2'],
    install_requires=[
        'appdirs',
        'bs4',
        'eventlet',
        'fiona',
        'haversine',
        'matplotlib',
        'netCDF4',
        'numpy',
        'ordered_set',
        'pandas',
        'paramiko',
        'psutil',
        # 'pygrib', # TODO: make separate install
        'pyproj>=2.6',
        'geopandas',
        'requests',
        'scipy',
        'seaborn',
        'shapely',
        'tropycal',
        'utm',
        'wget',
    ],
    # test and development dependencies
    extras_require={
        'testing': [
            'colored_traceback',
            'FileLock',
            'pytest',
            'pytest-mock',
            'pytest-cov',
            'pytest-xdist',
        ],
        'development': ['dunamai', 'flake8', 'isort', 'oitnb'],
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
