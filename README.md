# AdcircPy </h1>
## A Python interface for handling inputs and outputs for the ADCIRC hydrodynamic model. 

### Basic Installation:
You can install this software through conda or through pip.</br>
The dependecies of Python packages installed through pip must be satisfied by the operating system (for example, the headers for HDF5 for use with netCDF4). Therefore, the recommended method of installation is through conda. </br>

Method 1 describes the steps to install thorugh conda, while method 2 describes the steps to perform a pip only install.

#### Method 1: Conda (the easy way!)
Install Anaconda on your platform.</br>
Navigate to the directory where the AdcircPy has been downloaded and run the following commands:
```cmd
conda env create -f environment.yml
source activate AdcircPy
python setup.py install
```
This will install the AdcircPy on your conda environment. Remember to alwasy activate the AdcircPy enviroment by running the command
```cmd
source activate AdcircPy
```


#### Method 2: Pip-only (for the hardcore purist)

pyproj needs to be installed prior the installation on pip by running:</br>

```cmd
pip install git+https://github.com/jswhit/pyproj.git@master
```

Now you can install the package through pip by executing:
```cmd
pip install AdcircPy
```
If the package fails to install through pip, check if you have the headers installed in your os.

### Usage


This package is an API to handle input and output files for the [ADCIRC](http://adcirc.org) hydrodynamic model. 

Create a new direectory where you plan to start your project and create a new file with the following content:

```Python
from AdcircPy import AdcircPy
```

If it's the first time running AdircPy on this user account, the TPXO database will be cached to disk.

Now you can load an Adcirc mesh or an output file by doing: 

```Python
mesh = AdcircPy.read_mesh('/path/to/fort.14')
output = AdcircPy.read_output('/path/to/outputfile.nc'
```

The most basic you can do with this software is create plots:

```Python
import matplotlib.pyplot as plt 
mesh.make_plot()
output.make_plot()
plt.show()
```

Please report bugs to jreniel@gmail.com