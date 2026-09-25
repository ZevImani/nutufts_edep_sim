#!/bin/bash

# WE NEED TO SETUP ENV VARIABLES FOR ROOT, CUDA, OPENCV
alias python=python3
MACHINE=`uname --nodename`

# SETUP THE DEFAULT EDEPSIM DIRECTORY

# This is the install location of EdepSim, the core Geant4 code
export REPO_HOME_DIR=`pwd`/`dirname "${BASH_SOURCE[0]}"`/
export EDEPSIM_DIR=`pwd`/`dirname "${BASH_SOURCE[0]}"`/edep-sim/EDepSim
export EDEPSIM_BUILD_DIR=`pwd`/`dirname "${BASH_SOURCE[0]}"`/edep-sim/build
export EDEPSIM_SOURCE_DIR=`pwd`/`dirname "${BASH_SOURCE[0]}"`/edep-sim/edep-sim-source
export EDEPSIM_BIN_DIR=${EDEPSIM_DIR}/bin
export EDEPSIM_LIB_DIR=${EDEPSIM_DIR}/lib
export EDEPSIM_INCLUDE_DIR=${EDEPSIM_DIR}/include
export EDEPTLY_DIR=`pwd`/`dirname "${BASH_SOURCE[0]}"`/edeptly
export VECTORCLASS_INCLUDE_DIR=`pwd`/`dirname "${BASH_SOURCE[0]}"`/vectorclass
export SIMPLEDET_DIR=`pwd`/`dirname "${BASH_SOURCE[0]}"`/simpledet/
export SIMPLEDET_LIB_DIR=${SIMPLEDET_DIR}/build/lib/
export SIMPLEDET_PYTHON_DIR=${SIMPLEDET_DIR}/python/
export PETASTORM_DIR=`pwd`/`dirname "${BASH_SOURCE[0]}"`/petastorm/
export DARKNEWS_DIR=`pwd`/`dirname "${BASH_SOURCE[0]}"`/darknews/install/lib/python3.9/site-packages/

[[ ":$PATH:" != *":${EDEPSIM_BIN_DIR}:"* ]] && export PATH="${EDEPSIM_BIN_DIR}:${PATH}"
[[ ":$LD_LIBRARY_PATH:" != *":${EDEPSIM_LIB_DIR}:"* ]] && export LD_LIBRARY_PATH="${EDEPSIM_LIB_DIR}:${LD_LIBRARY_PATH}"
[[ ":$PYTHONPATH:" != *":${EDEPTLY_DIR}:"* ]] && export PYTHONPATH="${EDEPTLY_DIR}:${PYTHONPATH}"
[[ ":$PYTHONPATH:" != *":${SIMPLEDET_PYTHON_DIR}:"* ]] && export PYTHONPATH="${SIMPLEDET_PYTHON_DIR}:${PYTHONPATH}"
[[ ":$LD_LIBRARY_PATH:" != *":${SIMPLEDET_LIB_DIR}:"* ]] && export LD_LIBRARY_PATH="${SIMPLEDET_LIB_DIR}:${LD_LIBRARY_PATH}"
[[ ":$PYTHONPATH:" != *":${PETASTORM_DIR}:"* ]] && export PYTHONPATH="${PETASTORM_DIR}:${PYTHONPATH}"
[[ ":$PYTHONPATH:" != *":${DARKNEWS_DIR}:"* ]] && export PYTHONPATH="${DARKNEWS_DIR}:${PYTHONPATH}"

echo "SETUP FOR SINGULARITY CONTAINER"
source /usr/local/root/bin/thisroot.sh
source /usr/local/geant/geant4.10.06.p03/bin/geant4.sh

export CUDA_HOME=/usr/local/cuda/
[[ ":$LD_LIBRARY_PATH:" != *":${CUDA_HOME}/lib64:"* ]] && export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}"

export OPENCV_INCDIR=/usr/include
export OPENCV_LIBDIR=/usr/local/lib
    
echo "Key env variables set"
echo "REPO_HOM_DIR=${REPO_HOME_DIR}"
echo "EDEPSIM_DIR=${EDEPSIM_DIR}"
echo "SIMPLEDET_DIR=${SIMPLEDET_DIR}"
echo "PETASTORM_DIR=${PETASTORM_DIR}"
echo "DARKNEWS_DIR=${DARKNEWS_DIR}"

alias ls='ls --color=auto'
alias zlab='cd /cluster/tufts/wongjiradlabnu/zimani01'
alias zimani='cd /exp/uboone/app/users/imani'