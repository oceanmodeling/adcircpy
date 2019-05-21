# AdcircPy </h1>
### A Python interface for handling inputs and outputs for the ADCIRC hydrodynamic model. 

## Basic Installation:
You can install this software through conda or through pip.</br>
The dependecies of Python packages installed through pip must be satisfied by the operating system (for example, the headers for HDF5 for use with netCDF4). Therefore, the recommended method of installation is through conda. </br>

It is recommended that this package is used on Linux, but Windows should work as well. Not tested in MacOS, but should also work as well.</br>


Method 1 describes the steps to install thorugh conda, while method 2 describes the steps to perform a pip only install.

### Method 1: Conda (the easy way!)
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
#### Linux:
Same as Windows but use conda-linux-x86_64.yml instead.

### Method 2: Pip-only (for the hardcore purist)

pyproj needs to be installed prior the installation on pip by running:</br>

```cmd
pip install git+https://github.com/jswhit/pyproj.git@master
```

Now you can install the package through pip by executing:
```cmd
pip install AdcircPy
```
If the package fails to install through pip, check if you have the correct headers installed in your operating system. Using conda is highly recommended.

#### NOTE: If you are developer and would like to debug the package, you may install this package on developer mode by doing:
```cmd
pip install -e .
```

### Usage


This package is an API to handle input and output files for the [ADCIRC](http://adcirc.org) hydrodynamic model. 

Create a new directory where you plan to start your project and create a new file with the following content:

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
#### Tidal only run:

```Python
from AdcircPy import AdcircPy
from datetime import datetime, timedelta
Mesh = AdcircPy.read_mesh(fort14='/path/to/fort.14',
	                      fort13='/path/to/fort.13')
start_date = datetime.now()
end_date = start_date+timedelta(days=5)
tidalRun = Mesh.TidalRun(start_date, end_date)
tidalRun.dump('/directory/to/dump')
```

#### Generation of a Best Track Meteorological run:

```Python
from datetime import datetime, timedelta
from AdcircPy import AdcircPy
storm_id = 'AL182012'
spinup_date = datetime(2012, 10, 11, 0)
start_time  = spinup_date + timedelta(days=15)
end_time    = spinup_date + timedelta(days=19.25)
Mesh = AdcircPy.read_mesh(fort14='/path/to/fort.14',
	                      fort13='/path/to/fort.13')
BestTrackRun = Mesh.BestTrackRun('AL182012', start_time, end_time, spinup_date=spinup_date)
BestTrackRun.dump('/directory/to/dump')
```

#### Example where global outputs are requested:

```Python
from AdcircPy import AdcircPy
from AdcircPy import ElevationGlobalOutput as EGO
from datetime import datetime, timedelta
Mesh = AdcircPy.read_mesh(fort14='/path/to/fort.14',
	                      fort13='/path/to/fort.13')
start_date = datetime.now()
end_date = start_date+timedelta(days=5)
tidalRun = Mesh.TidalRun(start_date, end_date,
			ElevationGlobalOutput=EGO(sampling_frequency=timedelta(minutes=15)))
tidalRun.dump('/directory/to/dump')
```

See the "examples" directory for more examples where the Elevation stations outputs are requested.

Please report bugs to jreniel@gmail.com

pyproj, netCDF, scipy and gdal should be sourced from conda-forge. 