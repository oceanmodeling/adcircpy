#!/usr/bin/env python
import setuptools

conf = setuptools.config.read_configuration('./setup.cfg')
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
    # setup_requires=['pyproj', 'numpy'],
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
        'pyproj',
        'numpy',
        'ordered_set',
        # 'pygrib', # TODO: make separate install
        'psutil',
        'shapely'
    ],
    entry_points={
        'console_scripts': [

            # Generators
            'generate_hindcast=adcircpy.scripts.generate_hindcast:main',

            # Plotters
            'plot_mesh=adcircpy.scripts.plot_mesh:main',

            # 'PlotMaxele='
            # + 'adcircpy.scripts.PlotMaxele:main',

            # 'PlotElevationStationsOutput='
            # + 'adcircpy.scripts.PlotElevationStationsOutput:main',

            # 'PlotTidalStations='
            # + 'adcircpy.scripts.PlotElevationStationsOutput:main',

            'plot_fort61=adcircpy.scripts.plot_fort61:main',

            # Generators
            # 'GenerateTidalRun='
            # + 'adcircpy.scripts.GenerateTidalRun:main',

            'best_track_file=adcircpy.scripts.GenerateBestTrackFile:main',

            # 'GenerateBestTrackRun='
            # + 'adcircpy.scripts.GenerateBestTrackRun:main',

            # Validations
            # 'HighWaterMarkValidation='
            # + 'adcircpy.scripts.HighWaterMarkValidation:main',
        ]
    },
    test_suite='nose.collector',
    tests_require=['nose'])
