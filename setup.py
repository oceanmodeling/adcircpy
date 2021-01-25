#!/usr/bin/env python
import logging
import pathlib

import setuptools

try:
    try:
        from dunamai import Version
    except ImportError:
        import sys
        import subprocess

        subprocess.check_call(
                [sys.executable, '-m', 'pip', 'install', 'dunamai']
        )
        from dunamai import Version

    version = Version.from_any_vcs(
            pattern='^(?P<base>\d+\.\d+\.\d+)(-?((?P<stage>[a-zA-Z]+)\.?(?P<revision>\d+)?))?$'
    ).serialize()
except RuntimeError as error:
    logging.exception(error)
    logging.warning(f'using dummy version 0.0.0')
    version = '0.0.0'

metadata = setuptools.config.read_configuration(
        pathlib.Path(__file__).parent.absolute() / 'setup.cfg')['metadata']

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
            'development': [
                'colored_traceback',
                'coverage',
                'coverage-badge',
                'dunamai',
                'flake8',
                'nose',
                'rednose',
            ]
        },
        entry_points={
            'console_scripts': [
                # runtime commands
                'tidal_run=adcircpy.cmd.tidal_run:main',
                'best_track_run=adcircpy.cmd.best_track_run:main',
                # hack-a-mesh tools
                'iter_smooth=adcircpy.cmd.iter_smooth:main',
                # Generators
                # 'generate_hindcast=adcircpy._cmd.generate_hindcast:main',
                'best_track_file=adcircpy.cmd.best_track_file:main',
                # Plotters
                'plot_mesh=adcircpy.cmd.plot_mesh:main',
                'plot_maxele=adcircpy.cmd.plot_maxele:main',
                'plot_fort61=adcircpy.cmd.plot_fort61:main',
                'fort63=adcircpy.cmd.fort63:main',
                # tpxo_install
                'tpxo_install=adcircpy.forcing.tides.tpxo:install',
            ]
        },
        test_suite='nose.collector',
)
