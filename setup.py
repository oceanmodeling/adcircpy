import setuptools

# import os
# os.system('pip install cython')
# os.system('pip install git+https://github.com/jswhit/pyproj.git@master')

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="AdcircPy",
    version="0.8.0",
    author="Jaime R Calzada",
    author_email="jaime.calzada@gmail.com",
    description="Python package for working with ADCIRC input and output files.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jreniel/AdcircPy.git",
    # packages=['matplotlib',
    #                    'netCDF4',
    #                    'scipy',
    #                    'GDAL',
    #                    'haversine',
    #                    'cython',
    #                    'pyproj'],
    classifiers=[
        "Programming Language :: Python :: 3",
        # "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",],
    install_requires=['matplotlib',
                       'netCDF4',
                       'scipy',
                       'GDAL',
                       'haversine',
                       'cython',
                       'wget'],
    dependency_links=['git+https://github.com/jswhit/pyproj.git@master']

    # extras_require={"extras": [ "Cython==0.25.2", ],},
    # entry_points={'app_name':["foo = mypackage.some_module:foo",],}


)