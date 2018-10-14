

class PBS(object):
  def __init__(self, directory, **kwargs):
    self.directory = directory
    self._init_PBS_vars(**kwargs)

  def dump(self):
    with open(self.directory+'/qsub.pbs', 'w') as self.f:
      self.f.write('#!/bin/bash --login\n')
      self.f.write('#PBS -d .\n')
      self.__write_PBS_description()

  def _init_PBS_vars(**kwargs):
    

  def __write_PBS_description(self):
    if self.description is not None:
      f.write('#PBS -N {}'.format(self.description))



#PBS -A nosofs
#PBS -m be
#PBS -M jaime.calzadamarrero@upr.edu
#PBS -j oe
#PBS -o job.${PBS_JOBID}.oe
#PBS -l procs=600
#PBS -l walltime=01:45:00
#PBS -q batch

# Set path to ADCIRC binaries:
ADCIRC_BINARIES_PATH=/home/Jaime.Calzada/adcirc/bin/intel

# Set PATH to the location of the mesh file and fort.13:
FORT14_PATH=/home/Jaime.Calzada/work/HSOFS/CubaIkeModNOMAD1enoflux.grd
FORT13_PATH=/home/Jaime.Calzada/work/HSOFS/fort13nomad1elowerwaterdrag.grd

# Load HPC modules required to run ADCIRC
module load intel/15.6.233
module load impi/5.0.3.048
module load hdf5parallel
module load netcdf-hdf5parallel


#-------------Begin ADCIRC run---------------#
ROOTDIR=$PWD
WRITER_PROCS=1
RUN_PROCS="$(($PBS_NP-$WRITER_PROCS))"
COLDSTARTDIR=$ROOTDIR/coldstart
HOTSTARTDIR=$ROOTDIR/hotstart
OUTPUTDIR=$ROOTDIR/outputs

# Stage Coldstart
mkdir -p $COLDSTARTDIR
mkdir -p $HOTSTARTDIR
mkdir -p $OUTPUTDIR

ln -f $FORT14_PATH $COLDSTARTDIR/fort.14
ln -f $FORT13_PATH $COLDSTARTDIR/fort.13
ln -f $ROOTDIR/fort.15.coldstart $COLDSTARTDIR/fort.15

# Run coldstart
$ADCIRC_BINARIES_DIR/adcprep --np $RUN_PROCS --partmesh
$ADCIRC_BINARIES_DIR/adcprep --np $RUN_PROCS --prepall
mpirun -np $PBS_NP $ADCIRC_BINARIES_DIR/padcirc -W $WRITER_PROCS

# Stage hostart
ln -f $COLDSTARTDIR/fort.67.nc $HOTSTARTDIR/fort.67.nc
ln -f $COLDSTARTDIR/fort.15.hotstart $HOTSTARTDIR/fort.15
ln -f $FORT14_PATH $HOTSTARTDIR/fort.14
ln -f $FORT13_PATH $HOTSTARTDIR/fort.13

# Prepare Best Track wind file.
# cp $ROOTDIR/fort.22.best_track $ROOTDIR/fort.22
# $ADCIRC_BINARIES_DIR/aswip -n 20 -m 4 -z 2
# rm $ROOTDIR/fort.22
# ln -f ./NWS_20_fort.22 ./fort.22

# run hotstart
$ADCIRC_BINARIES_DIR/adcprep --np $RUN_PROCS --partmesh
$ADCIRC_BINARIES_DIR/adcprep --np $RUN_PROCS --prepall
mpirun -np $PBS_NP $ADCIRC_BINARIES_DIR/padcirc -W $WRITER_PROCS

# cleanup
ln -f $HOTSTARTDIR/*.nc $OUTPUTDIR/
