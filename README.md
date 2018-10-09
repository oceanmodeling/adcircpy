[![CircleCI](https://circleci.com/gh/jreniel/AdcircPy/tree/master.svg?style=svg)](https://circleci.com/gh/jreniel/AdcircPy/tree/master)

# AdcircPy </h1>
## A Python interface for handling inputs and outputs for the ADCIRC hydrodynamic model. 

### Basic Installation:
You can install this software through conda or through pip.</br>
The dependecies of Python packages installed through pip must be satisfied by the operating system (for example, the headers for HDF5 for use with netCDF4). Therefore, the recommended method of installation is through conda. </br>

Method 1 describes the steps to install thorugh conda, while method 2 describes the steps to perform a pip only install.

#### Method 1: Conda (the easy way!)
#### Windows:
Install Anaconda for Windows</br>
Using the Anaconda Command Prompt, navigate to the directory where the AdcircPy has been downloaded and run the following command:
```cmd
conda env create -f conda-windows-x86_64.yml
```
After the installation of the dependecies is done, you should activate the AdcircPy environment by running:
```cmd
conda activate AdcircPy
```
Finally, you can install the package by issuing the command:
```cmd
python setup.py install
```
This will install the AdcircPy on your conda environment. Remember to alwasy activate the AdcircPy enviroment by running the command
```cmd
conda activate AdcircPy
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
output = AdcircPy.read_output('/path/to/outputfile.nc')
```

The most basic you can do with this software is create plots:

```Python
import matplotlib.pyplot as plt 
mesh.make_plot()
output.make_plot()
plt.show()
```

### Example of fort.15 generation:
This is still in development, so it will only output a partial fort.15
```Python
from AdcircPy import AdcircPy
from AdcircPy import TidalForcing
from AdcircPy import AdcircRun
from datetime import datetime, timedelta

Mesh = AdcircPy.read_mesh(fort14='/path/to/fort.14',
	                      fort13='/path/to/fort.13')
start_date = datetime.now()
end_date = start_date+timedelta(days=5)
tidalForcing = TidalForcing(start_date, end_date,
	# optionally may pass a constituent list for forcing
							constituents=None, 
	# Optionally pass a spinup_date
	# The package will take 15 days prior to start_date as spinup by default.
							spinup_date=None)
adcircRun = AdcircRun(Mesh, Tides=tidalForcing)
# Writes the files to directory
adcircRun.dump('/directory/to/dump')
```


Please report bugs to jreniel@gmail.com