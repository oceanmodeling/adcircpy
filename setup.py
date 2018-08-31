import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="AdcircPy",
    version="0.8.0",
    author="Jaime R Calzada",
    author_email="jaime.calzada@noaa.gov",
    description="Python package for working with ADCIRC input and output files.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jreniel/AdcircPy.git",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        # "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)