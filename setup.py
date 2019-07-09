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
            'PlotMesh=adcircpy.scripts.PlotMesh:main',
            'GenerateTidalRun=adcircpy.scripts.GenerateTidalRun:main',
            'PlotTidalStationsOutput=adcircpy.scripts.PlotTidalStationsOutput:main',  # noqa:E501
            'GenerateBestTrackFile=adcircpy.scripts.GenerateBestTrackFile:main',  # noqa:E501
            'PlotMaxele=adcircpy.scripts.PlotMaxele:main',
            'GenerateBestTrackRun=adcircpy.scripts.GenerateBestTrackRun:main',  # noqa:E501
            'HighWaterMarkValidation=adcircpy.scripts.HighWaterMarkValidation:main',  # noqa:E501
        ]
    },
    test_suite='nose.collector',
    tests_require=['nose'])
