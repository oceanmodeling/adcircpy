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
    setup_requires=['wheel'],
    install_requires=[
        'numpy',
        'matplotlib',
        'netCDF4',
        'scipy',
        'haversine',
        'wget',
        'utm',
        'gdal',
        'requests',
        'tqdm',
        'eventlet',
        'bs4',
        'seaborn',
        'pandas',
        'pyproj',
        'pygrib'
    ],
    entry_points={
        'console_scripts': [

            # Plotters
            'PlotMesh='
            + 'adcircpy.scripts.PlotMesh:main',

            'PlotMaxele='
            + 'adcircpy.scripts.PlotMaxele:main',

            'PlotElevationStationsOutput='
            + 'adcircpy.scripts.PlotElevationStationsOutput:main',

            'PlotTidalStations='
            + 'adcircpy.scripts.PlotElevationStationsOutput:main',

            'PlotFort61='
            + 'adcircpy.scripts.PlotElevationStationsOutput:main',

            # Generators
            'GenerateTidalRun='
            + 'adcircpy.scripts.GenerateTidalRun:main',

            'GenerateBestTrackFile='
            + 'adcircpy.scripts.GenerateBestTrackFile:main',

            'GenerateBestTrackRun='
            + 'adcircpy.scripts.GenerateBestTrackRun:main',

            # Validations
            'HighWaterMarkValidation='
            + 'adcircpy.scripts.HighWaterMarkValidation:main',
        ]
    },
    test_suite='nose.collector',
    tests_require=['nose'])
