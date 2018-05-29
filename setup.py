import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="AdcircPy",
    version="0.0.0",
    author="Jaime R Calzada",
    author_email="jreniel@dalidev.com",
    description="Python tools for working with the ADCIRC model.",
    # long_description=long_description,
    # long_description_content_type="text/markdown",
    # url="https://github.com/pypa/example-project",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        # "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
)