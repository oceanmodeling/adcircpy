#!/usr/bin/env python
import json
import pathlib

import setuptools

try:
    import requests
except ImportError:
    import sys
    import subprocess

    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'requests'])
    import requests

try:
    from dunamai import Version
except ImportError:
    import sys
    import subprocess

    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'dunamai'])
    from dunamai import Version

try:
    from packaging.version import parse as package_parse
except ImportError:
    from pip._vendor.packaging.version import parse as package_parse


def latest_pypi_version(package: str, prerelease: bool = False) -> str:
    """ https://stackoverflow.com/a/34366589 """
    response = requests.get(f'https://pypi.python.org/pypi/{package}/json')

    if response.status_code == requests.codes.ok:
        releases = json.loads(response.text).get('releases', [])
        versions = list(sorted(package_parse(release) for release in releases))
        if not prerelease:
            versions = [version for version in versions
                        if not version.is_prerelease]
        version = str(max(versions))
    else:
        version = None
    return version


PARENT = pathlib.Path(__file__).parent.absolute()
conf = setuptools.config.read_configuration(PARENT / 'setup.cfg')
meta = conf['metadata']

try:
    version = Version.from_any_vcs(
            pattern='^(?P<base>\d+\.\d+\.\d+)(-?((?P<stage>[a-zA-Z]+)\.?(?P<revision>\d+)?))?$'
    ).serialize()
except:
    # TODO: this only retrieves the latest version from PyPI; the user might have requested a different version
    version = latest_pypi_version('adcircpy')

setuptools.setup(
        name=meta['name'],
        version=version,
        author=meta['author'],
        author_email=meta['author_email'],
        description=meta['description'],
        long_description=meta['long_description'],
        long_description_content_type='text/markdown',
        url=meta['url'],
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
