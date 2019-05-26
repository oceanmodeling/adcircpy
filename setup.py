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
        'pyproj'
    ],
    entry_points={
        'console_scripts': [
            'PlotMesh=AdcircPy.entrypoints.PlotMesh:main',
            'GenerateTidalRun=AdcircPy.entrypoints.GenerateTidalRun:main',
            'PlotTidalStationsOutput=AdcircPy.entrypoints.PlotTidalStationsOutput:main',  # noqa:E501
            'GenerateBestTrackFile=AdcircPy.entrypoints.GenerateBestTrackFile:main',  # noqa:E501
            'PlotMaxele=AdcircPy.entrypoints.PlotMaxele:main',
            'GenerateBestTrackRun=AdcircPy.entrypoints.GenerateBestTrackRun:main',  # noqa:E501
            'HighWaterMarkValidation=AdcircPy.entrypoints.HighWaterMarkValidation:main',  # noqa:E501
        ]
    },
    test_suite='nose.collector',
    tests_require=['nose'])
