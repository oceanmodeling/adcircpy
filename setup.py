#!/usr/bin/env python
import setuptools
import pathlib
module_path = pathlib.Path(__file__).parent.absolute()
conf = setuptools.config.read_configuration(module_path / 'setup.cfg')
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
            'generate_hindcast=adcircpy.cmd.generate_hindcast:main',
            'best_track_file=adcircpy.cmd.GenerateBestTrackFile:main',

            # Plotters
            'plot_mesh=adcircpy.cmd.plot_mesh:main',
            'plot_maxele=adcircpy.cmd.plot_maxele:main',
            'plot_fort61=adcircpy.cmd.plot_fort61:main',

        ]
    },
    test_suite='nose.collector',
    tests_require=['nose'])
