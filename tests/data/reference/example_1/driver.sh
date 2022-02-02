#!/bin/bash

ulimit -s unlimited

set -e

NPROCS=4

main() {
  rm -rf work
  mkdir work
  cd work
  ln -sf ../fort.14
  ln -sf ../fort.13
  ln -sf ../fort.15 ./fort.15
  adcprep --np 4 --partmesh
  adcprep --np 4 --prepall
  mpiexec -n 4 padcirc
  clean_directory
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
