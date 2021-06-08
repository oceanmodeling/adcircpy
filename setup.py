#!/usr/bin/env python
import importlib
import logging
import os
from pathlib import Path
import subprocess
import sys

import setuptools

BUILT_PACKAGES = {
    'fiona': ['gdal'],
    'matplotlib': [],
    'netCDF4': [],
    'numpy': [],
    'pyproj': ['gdal'],
    'scipy': [],
    'shapely': ['gdal'],
}
is_conda = (Path(sys.prefix) / 'conda-meta').exists()

if is_conda:
    conda_packages = []
    for conda_package in BUILT_PACKAGES:
        try:
            importlib.import_module(conda_package)
        except:
            conda_packages.append(conda_package)
    if len(conda_packages) > 0:
        subprocess.check_call(['conda', 'install', '-y', *conda_packages])

if os.name == 'nt':
    for required_package, pipwin_dependencies in BUILT_PACKAGES.items():
        try:
            importlib.import_module(required_package)
        except:
            try:
                import pipwin
            except:
                subprocess.check_call(
                        [sys.executable, '-m', 'pip', 'install', 'pipwin'])

            failed_pipwin_packages = []
            for pipwin_package in pipwin_dependencies + [required_package]:
                try:
                    subprocess.check_call(
                            [sys.executable, '-m', 'pipwin', 'install',
                             pipwin_package.lower()]
                    )
                except subprocess.CalledProcessError:
                    failed_pipwin_packages.append(pipwin_package)

            if len(failed_pipwin_packages) > 0:
                raise RuntimeError(
                        f'failed to download or install non-conda Windows build(s) of {" and ".join(failed_pipwin_packages)}; you can either\n'
                        '1) install within an Anaconda environment, or\n'
                        f'2) `pip install <file>.whl`, with `<file>.whl` downloaded from {" and ".join("https://www.lfd.uci.edu/~gohlke/pythonlibs/#" + value.lower() for value in failed_pipwin_packages)} for your Python version'
                )

try:
    try:
        from dunamai import Version
    except ImportError:
        import subprocess
        import sys

        subprocess.check_call(
                [sys.executable, '-m', 'pip', 'install', 'dunamai'])
        from dunamai import Version

    version = Version.from_any_vcs().serialize()
except RuntimeError as error:
    logging.exception(error)
    version = '0.0.0'

logging.info(f'using version {version}')

metadata = setuptools.config.read_configuration(
        Path(__file__).parent.absolute() / 'setup.cfg'
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
