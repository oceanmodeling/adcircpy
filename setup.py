from setuptools import setup, find_packages
from setuptools.config import read_configuration
conf = read_configuration('./setup.cfg')
meta=conf['metadata']
setup(name=meta['name'],
    version=meta['version'],
    author=meta['author'],
    author_email=meta['author_email'],
    description=meta['description'],
    long_description=meta['long_description'],
    long_description_content_type="text/markdown",
    url=meta['url'],
    packages=find_packages(),
    setup_requires=['cython'],
    install_requires=[ 'numpy',
                       'matplotlib',
                       'netCDF4',
                       'scipy',
                       'GDAL',
                       'haversine',
                       'pyproj',
                       'wget',
                       'utm',
                       'requests'],
    test_suite='nose.collector',
    tests_require=['nose'])