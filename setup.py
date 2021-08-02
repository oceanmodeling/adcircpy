import logging
import os
from pathlib import Path
import re
import subprocess
import sys

from setuptools import config, find_packages, setup

DEPENDENCIES = {
    'appdirs': [],
    'bs4': [],
    'eventlet': [],
    'fiona': ['gdal'],
    'haversine': [],
    'matplotlib': [],
    'netCDF4': [],
    'numpy': [],
    'ordered_set': [],
    'pandas': [],
    'paramiko': [],
    'psutil': [],
    # 'pygrib': [], # TODO: make separate install
    'pyproj>=2.6': [],
    'geopandas': [],
    'requests': [],
    'scipy': [],
    'seaborn': [],
    'shapely': [],
    'tropycal': [],
    'utm': [],
    'wget': [],
}


def installed_packages() -> [str]:
    return [
        re.split('#egg=', re.split('==| @ ', package.decode())[0])[-1].lower()
        for package in subprocess.check_output(
                [sys.executable, '-m', 'pip', 'freeze']
        ).splitlines()
    ]


def missing_packages(dependencies: {str: [str]}) -> {str: [str]}:
    return {
        dependency: subdependencies
        for dependency, subdependencies in dependencies.items()
        if re.split('<|<=|==|>=|>', dependency)[
               0].lower() not in installed_packages()
    }


missing_dependencies = missing_packages(DEPENDENCIES)

if (Path(sys.prefix) / 'conda-meta').exists() and len(
        missing_dependencies) > 0:
    try:
        subprocess.check_call(
                ['conda', 'install', '-y', list(missing_dependencies)])
    except:
        for dependency in list(missing_dependencies):
            try:
                subprocess.check_call(['conda', 'install', '-y', dependency])
            except:
                continue

    missing_dependencies = missing_packages(DEPENDENCIES)

if os.name == 'nt' and len(missing_dependencies) > 0:
    if 'pipwin' not in installed_packages():
        subprocess.check_call(
                [sys.executable, '-m', 'pip', 'install', 'pipwin'])
    subprocess.check_call([sys.executable, '-m', 'pipwin', 'refresh'])

    for dependency, subdependencies in missing_dependencies.items():
        failed_pipwin_packages = []
        for _ in range(1 + len(subdependencies)):
            for package_name in [dependency] + subdependencies:
                if package_name in missing_packages(DEPENDENCIES):
                    try:
                        subprocess.check_call(
                                [sys.executable, '-m', 'pipwin', 'install',
                                 package_name.lower()]
                        )
                        if package_name in failed_pipwin_packages:
                            failed_pipwin_packages.remove(package_name)
                    except subprocess.CalledProcessError:
                        failed_pipwin_packages.append(package_name)

            # since we don't know the dependencies here, repeat this process n number of times
            # (worst case is `O(n)`, where the first package is dependant on all the others)
            if len(failed_pipwin_packages) == 0:
                break

    missing_dependencies = missing_packages(DEPENDENCIES)

try:
    if 'dunamai' not in installed_packages():
        subprocess.check_call(
                [sys.executable, '-m', 'pip', 'install', 'dunamai'])

    from dunamai import Version

    version = Version.from_any_vcs().serialize()
except RuntimeError as error:
    logging.exception(error)
    version = '0.0.0'

logging.info(f'using version {version}')

metadata = config.read_configuration('setup.cfg')['metadata']

setup(
        name=metadata['name'],
        version=version,
        author=metadata['author'],
        author_email=metadata['author_email'],
        description=metadata['description'],
        long_description=metadata['long_description'],
        long_description_content_type='text/markdown',
        url=metadata['url'],
        packages=find_packages(),
        python_requires='>=3.6',
        setup_requires=['dunamai', 'requests', 'setuptools>=41.2'],
        install_requires=list(DEPENDENCIES),
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
