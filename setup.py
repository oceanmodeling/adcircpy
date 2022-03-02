import warnings

from dunamai import Version
from setuptools import find_packages, setup

try:
    __version__ = Version.from_any_vcs().serialize()
except RuntimeError as error:
    warnings.warn(f'{error.__class__.__name__} - {error}')
    __version__ = '0.0.0'

setup(version=__version__, packages=find_packages(exclude=('tests',)), test_suite='tests')
