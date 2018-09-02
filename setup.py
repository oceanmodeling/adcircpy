from setuptools import setup
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
    setup_requires=['cython', 'numpy'],
    install_requires=[ 'matplotlib',
                       'netCDF4',
                       'scipy',
                       'GDAL',
                       'haversine',
                       'pyproj',
                       'wget']
)