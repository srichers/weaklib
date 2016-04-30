#!/bin/bash

export WEAKLIB_MACHINE=$1

if [[ $WEAKLIB_MACHINE == titan* ]]; then

  echo
  echo "INFO: Setting environment for" $WEAKLIB_MACHINE

  source ${MODULESHOME}/init/bash

  module unload fftw cray-hdf5 cray-petsc silo subversion
  module unload pgi gcc cce pathscale
  module unload PrgEnv-pgi PrgEnv-gnu PrgEnv-cray PrgEnv-pathscale PrgEnv-intel

elif [[ $WEAKLIB_MACHINE == sn1987b* ]]; then

  echo
  echo "INFO: Setting environment for" $WEAKLIB_MACHINE

fi


if [[ $WEAKLIB_MACHINE == titan_gnu ]]; then

  echo

  module load PrgEnv-gnu
  module load cray-hdf5

elif [[ $WEAKLIB_MACHINE == sn1987b ]]; then

  echo

else

  echo "  WARNING: Unknown machine " $WEAKLIB_MACHINE

fi