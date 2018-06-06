This package was written with Python 3 in mind. It can potentially work in Python 2, but there are no guarantees.

To install dependencies: <br \> 

This package uses GDAL for Python (Geospatial Data Abstraction Library). <br \>
The easiest way to get GDAL installed is through conda (preferrably Minicoda). <br \>
Installing Miniconda has the added advantage of a neat Python environment manager with a minimal footprint.
Once conda is installed, you can do
```bash
conda env create -f adcpy.yml
```

This will create a new Python environment named adcpy, with all the dependencies included.<br \>
You can activate this environment with:
```bash
source activate adcpy
```

Since this package has not been officially publised to the conda or pip repositories you will need to install it manually, but this is 
very easy to do! <br />

Just cd inside the directory you just cloned (named AdcircPythonTools by default) and do:
```bash
echo "PYTHONPATH=$PWD:$PYTHONPATH" >> ~/.bashrc && echo "PATH=$PWD/apps:$PATH" >> ~/.bashrc
```

Now you can run all the standalone apps in the apps folder by just invoking their name in the command line.<br \>

The apps include a help command you can do 
```bash
plotMesh -h
```

You can also import the code in your Python scripts by doing
```python
import AdcircPy
```

Note that this is an unfinished product, and while most of it will work, it is for demostration purposes only!
Most of the standalone apps are unfinished but the AdcircPy module is in an advanced stage.

For questions: <br />
jreniel@gmail.com


