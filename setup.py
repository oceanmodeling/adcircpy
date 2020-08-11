#!/usr/bin/env python
import pathlib

import setuptools

PARENT = pathlib.Path(__file__).parent.absolute()
conf = setuptools.config.read_configuration(PARENT / 'setup.cfg')
meta = conf['metadata']
setuptools.setup(
    name=meta['name'],
    version=meta['version'],
    author=meta['author'],
    author_email=meta['author_email'],
    description=meta['description'],
    long_description=meta['long_description'],
    long_description_content_type="text/markdown",
    url=meta['url'],
    packages=setuptools.find_packages(),
    python_requires='>=3.6',
    setup_requires=["setuptools>=41.2"],
    install_requires=[
        'matplotlib',
        'netCDF4',
        'scipy',
        'haversine',
        'wget',
        'utm',
        'requests',
        'eventlet',
        'bs4',
        'seaborn',
        'pandas',
        'pyproj>=2.6',
        'numpy',
        'ordered_set',
        # 'pygrib', # TODO: make separate install
        'psutil',
        'shapely',
        'paramiko',
        'fiona',
        'wget',
        'tropycal',
        # test and development dependencies
        'nose',
        'rednose',
        'coverage-badge',
        'colored_traceback',
    ],
    entry_points={
        'console_scripts': [

            # runtime commands
            "tidal_run=adcircpy.cmd.tidal_run:main",
            "best_track_run=adcircpy.cmd.best_track_run:main",

            # hack-a-mesh tools
            "iter_smooth=adcircpy.cmd.iter_smooth:main",

            # Generators
            # 'generate_hindcast=adcircpy._cmd.generate_hindcast:main',
            'best_track_file=adcircpy.cmd.best_track_file:main',

            # Plotters
            'plot_mesh=adcircpy.cmd.plot_mesh:main',
            'plot_maxele=adcircpy.cmd.plot_maxele:main',
            'plot_fort61=adcircpy.cmd.plot_fort61:main',
            'fort63=adcircpy.cmd.fort63:main',

            # tpxo_install
            'tpxo_install=adcircpy.forcing.tides.tpxo:install'

        ]
    },
    test_suite='nose.collector',
    tests_require=['nose'])
