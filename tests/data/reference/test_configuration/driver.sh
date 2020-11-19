#!/bin/bash

ulimit -s unlimited

set -e

NPROCS=4

main() {
  SECONDS=0
  run_coldstart_phase
  if grep -Rq "ERROR: Elevation.gt.ErrorElev, ADCIRC stopping." padcirc.log; then
    duration=$SECONDS
    echo "ERROR: Elevation.gt.ErrorElev, ADCIRC stopping."
    echo "Wallclock time: $(($duration / 60)) minutes and $(($duration % 60)) seconds."
    exit -1
  else
    run_hotstart_phase
    duration=$SECONDS
    if grep -Rq "ERROR: Elevation.gt.ErrorElev, ADCIRC stopping." padcirc.log; then
      echo "ERROR: Elevation.gt.ErrorElev, ADCIRC stopping."
      echo "Wallclock time: $(($duration / 60)) minutes and $(($duration % 60)) seconds."
      exit -1
    fi
  fi
  echo "Wallclock time: $(($duration / 60)) minutes and $(($duration % 60)) seconds."
}

run_coldstart_phase() {
  rm -rf coldstart
  mkdir coldstart
  cd coldstart
  ln -sf ../fort.14
  ln -sf ../fort.13
  ln -sf ../fort.15.coldstart ./fort.15
  adcprep --np 4 --partmesh
  adcprep --np 4 --prepall
  mpiexec -n 4 padcirc 2>&1 | tee ../padcirc.logclean_directory
  cd ..
}

run_hotstart_phase() {
  rm -rf hotstart
  mkdir hotstart
  cd hotstart
  ln -sf ../fort.14
  ln -sf ../fort.13
  ln -sf ../fort.15.hotstart ./fort.15
  ln -sf ../coldstart/fort.67.nc
  ln -sf ../fort.22 ./fort.22
  aswip
  mv NWS_17_fort.22 fort.22
  adcprep --np 4 --partmesh
  adcprep --np 4 --prepall
  mpiexec -n 4 padcirc 2>&1 | tee -a ../padcirc.logclean_directory
  cd ..
}

clean_directory() {
  rm -rf PE*
  rm -rf partmesh.txt
  rm -rf metis_graph.txt
  rm -rf fort.13
  rm -rf fort.14
  rm -rf fort.15
  rm -rf fort.16
  rm -rf fort.80
  rm -rf fort.68.nc
}

main
