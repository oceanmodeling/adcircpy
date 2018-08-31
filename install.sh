#!/bin/bash -e
############################
# Setup script for AdcircPy#
############################

# while [[ "$#" > 0 ]]; do case $1 in
#   -d|--deploy) deploy="$2"; shift;;
#   -u|--uglify) uglify=1;;
#   *) echo "Unknown parameter passed: $1"; exit 1;;
# esac; shift; done

GIT_DIR=$(dirname `which $0`)
cd $GIT_DIR

# Installation Directory
INSTALLDIR=~/.local/Hurricaned

# environment file
RCFILE=$INSTALLDIR/.rcfile
BASH_ALIAS=$INSTALLDIR/.bash_alias

# Check if wget exists.
if [ ! -x /usr/bin/wget ] ; then
    command -v wget >/dev/null 2>&1 || { echo >&2 "Please install wget or set it in your path. Aborting."; exit 1; }
fi

# Create installation directory
mkdir -p $INSTALLDIR

# Install conda dependencies
if [ ! -d $INSTALLDIR/miniconda3 ]; then
	wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -P $INSTALLDIR
	bash $INSTALLDIR/Miniconda3-latest-Linux-x86_64.sh -b -p $INSTALLDIR/miniconda3
	source  $INSTALLDIR/miniconda3/bin/activate $INSTALLDIR/miniconda3
	conda create -y --name Hurricaned python=3
	conda update -y -n base conda
	conda activate Hurricaned
	pip install --upgrade pip
	pip install matplotlib scipy haversine requests bs4 psycopg2-binary
	conda install -y -c anaconda netcdf4 pyproj gdal
	if grep -Fxq "source $BASH_ALIAS" ~/.bashrc
		then
		    true
		else
			echo "# Line added by Hurricaned installer" >> ~/.bashrc
			echo "source $BASH_ALIAS" >> ~/.bashrc
	fi
	if grep -Fxq "source $BASH_ALIAS" ~/.zshrc
	then
	    true
	else
		echo "# Line added by Hurricaned installer" >> ~/.zshrc
		echo "source $BASH_ALIAS" >> ~/.zshrc
	fi
fi

rsync -za --delete-before $GIT_DIR/AdcircPy/ $INSTALLDIR/AdcircPy/
rsync -za --delete-before $GIT_DIR/bin/ $INSTALLDIR/bin/

if [ ! -f $INSTALLDIR/AdcircPy/core/h_tpxo9.v1.nc ]; then
	if [ ! -f $INSTALLDIR/tpxo9_netcdf.tar.gz ]; then
		wget -c ftp://ftp.oce.orst.edu/dist/tides/Global/tpxo9_netcdf.tar.gz -P $INSTALLDIR
	fi
	echo "Extracting TPXO file..."
	tar -xf $INSTALLDIR/tpxo9_netcdf.tar.gz -C $INSTALLDIR/AdcircPy/core h_tpxo9.v1.nc
fi

echo "
source $INSTALLDIR/miniconda3/bin/activate $INSTALLDIR/miniconda3
source activate Hurricaned
export PATH=$INSTALLDIR/bin:\$PATH
export PYTHONPATH="${PYTHONPATH}:$INSTALLDIR"
" > $RCFILE
echo "alias Hurricaned=\"source $RCFILE\"" > $BASH_ALIAS
echo "alias hurricaned.list='ls $INSTALLDIR/bin'" >> $BASH_ALIAS

echo "
##################################################################
# PLEASE READ THE FOLLOWING NOTICE:                              #
##################################################################

The installation of Hurricaned was succesfull!
Activate the environment by running the command \"Hurricaned\"

$ Hurricaned

If this is your first time installing, and Hurricaned says \"command not found\",
try resourcing your ~/.bashrc file first by running the following command:

$ source ~/.bashrc && Hurricaned

Happy Modeling!
"