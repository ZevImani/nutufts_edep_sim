#!/bin/bash

# WE NEED TO SETUP ENV VARIABLES FOR ROOT, CUDA, OPENCV
alias python=python3
MACHINE=`uname --nodename`

# SETUP THE DEFAULT EDEPSIM DIRECTORY

# This is the install location of EdepSim, the core Geant4 code.
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

[[ ":$PATH:" != *":${EDEPSIM_BIN_DIR}:"* ]] && export PATH="${EDEPSIM_BIN_DIR}:${PATH}"
[[ ":$LD_LIBRARY_PATH:" != *":${EDEPSIM_LIB_DIR}:"* ]] && export LD_LIBRARY_PATH="${EDEPSIM_LIB_DIR}:${LD_LIBRARY_PATH}"
[[ ":$PYTHONPATH:" != *":${EDEPTLY_DIR}:"* ]] && export PYTHONPATH="${EDEPTLY_DIR}:${PYTHONPATH}"
[[ ":$PYTHONPATH:" != *":${SIMPLEDET_PYTHON_DIR}:"* ]] && export PYTHONPATH="${SIMPLEDET_PYTHON_DIR}:${PYTHONPATH}"
[[ ":$LD_LIBRARY_PATH:" != *":${SIMPLEDET_LIB_DIR}:"* ]] && export LD_LIBRARY_PATH="${SIMPLEDET_LIB_DIR}:${LD_LIBRARY_PATH}"
[[ ":$PYTHONPATH:" != *":${PETASTORM_DIR}:"* ]] && export PYTHONPATH="${PETASTORM_DIR}:${PYTHONPATH}"

if [ $MACHINE == "trex" ]
then
    echo "SETUP TREX"
    source /usr/local/root/6.16.00_py3/bin/thisroot.sh

    export CUDA_HOME=/usr/local/cuda-10.0
    [[ ":$LD_LIBRARY_PATH:" != *":${CUDA_HOME}/lib64:"* ]] && export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}"

    export OPENCV_INCDIR=/usr/local/opencv/opencv-3.4.6/include
    export OPENCV_LIBDIR=/usr/local/opencv/opencv-3.4.6/lib

elif [ $MACHINE == "meitner" ]
then
    echo "SETUP MEITNER"

    source /usr/local/bin/thisroot.sh

    export CUDA_HOME=/usr/local/cuda-10.0
    [[ ":$LD_LIBRARY_PATH:" != *":${CUDA_HOME}/lib64:"* ]] && export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}"

    export OPENCV_INCDIR=/usr/include/opencv4
    export OPENCV_LIBDIR=/usr/lib/x86_64-linux-gnu

elif [ $MACHINE == "goeppert" ]
then
    echo "SETUP GOEPPERT"
    source /usr/local/root/6.16.00_py3/bin/thisroot.sh

    export CUDA_HOME=/usr/local/cuda/
    [[ ":$LD_LIBRARY_PATH:" != *":${CUDA_HOME}/lib64:"* ]] && export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}"

    export OPENCV_INCDIR=/usr/local/include
    export OPENCV_LIBDIR=/usr/local/lib    

elif [ $MACHINE == "mayer" ]
then
    echo "SETUP MAYER"
    source /usr/local/root/6.16.00_python3/bin/thisroot.sh

    export CUDA_HOME=/usr/local/cuda-10.0
    [[ ":$LD_LIBRARY_PATH:" != *":${CUDA_HOME}/lib64:"* ]] && export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}"

    export OPENCV_INCDIR=/usr/local/include
    export OPENCV_LIBDIR=/usr/local/lib    

elif [ $MACHINE == "blade" ]
then
    echo "SETUP TARITREE's RAZER BLADE"
    source /home/twongjirad/software/root6/py3_build/bin/thisroot.sh

    export CUDA_HOME=/usr/local/cuda
    [[ ":$LD_LIBRARY_PATH:" != *":${CUDA_HOME}/lib64:"* ]] && export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}"

    export OPENCV_INCDIR=/usr/include
    export OPENCV_LIBDIR=/usr/local/lib    
    
elif [ $MACHINE == "pop-os" ]
then
    echo "SETUP TARITREE's RAZER BLADE PRO"
    #source /home/twongjirad/software/root/build_gcc11/bin/thisroot.sh
    #source /usr/local/root_6.32.02/bin/thisroot.sh
    source /usr/local/root_v6.28.12_py3.10/bin/thisroot.sh

    export CUDA_HOME=/usr/local/cuda-11.6
    [[ ":$LD_LIBRARY_PATH:" != *":${CUDA_HOME}/lib64:"* ]] && export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}"
    [[ ":$PATH:" != *":${CUDA_HOME}/bin:"* ]] && export PATH="${CUDA_HOME}/bin:${PATH}"

    export OPENCV_INCDIR=/usr/include/opencv4
    export OPENCV_LIBDIR=/usr/lib/x86_64-linux-gnu/

    # setup geant4
    source /usr/local/geant4/g4_10.6.3/bin/geant4.sh
    
elif [ $MACHINE == "mmr-Alienware-x15-R1" ]
then
    echo "SETUP MATT's ALIENWARE x15 R1"
    source /home/matthew/software/root/install/bin/thisroot.sh

    export CUDA_HOME=/usr/local/cuda-11.3/
    [[ ":$LD_LIBRARY_PATH:" != *":${CUDA_HOME}/lib64:"* ]] && export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}"

    export OPENCV_INCDIR=/usr/include/opencv4
    export OPENCV_LIBDIR=/usr/lib/x86_64-linux-gnu/

    export LIBTORCH_DIR="/home/matthew/software/pythonVEnvs/myPy3.8.10/lib/python3.8/site-packages/torch"
    export LIBTORCH_LIBDIR=${LIBTORCH_DIR}/lib
    export LIBTORCH_INCDIR=${LIBTORCH_DIR}/include
    [[ ":$LD_LIBRARY_PATH:" != *":${LIBTORCH_LIBDIR}:"* ]] && \
        export LD_LIBRARY_PATH="${LIBTORCH_LIBDIR}:${LD_LIBRARY_PATH}"
   
else
    echo "DEFAULT SETUP (COMPAT WITH SINGULARITY CONTAINER)"
    source /usr/local/root/bin/thisroot.sh
    source /usr/local/geant/geant4.10.06.p03/bin/geant4.sh
    export JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64

    export CUDA_HOME=/usr/local/cuda/
    [[ ":$LD_LIBRARY_PATH:" != *":${CUDA_HOME}/lib64:"* ]] && export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}"

    export OPENCV_INCDIR=/usr/include
    export OPENCV_LIBDIR=/usr/local/lib
    
fi


# LIBTORCH
# location below is typically where running `pip install torch` will put pytorch
# export LIBTORCH_DIR="/usr/local/lib/python2.7/dist-packages/torch"
#export LIBTORCH_DIR="/home/jmills/.local/lib/python3.5/site-packages/torch"
#export LIBTORCH_LIBDIR=${LIBTORCH_DIR}/lib
#export LIBTORCH_INCDIR=${LIBTORCH_DIR}/lib/include
#export LIBTORCH_DIR="/usr/local/torchlib"
#export LIBTORCH_LIBDIR=${LIBTORCH_DIR}/lib
#export LIBTORCH_INCDIR=${LIBTORCH_DIR}/include
#[[ ":$LD_LIBRARY_PATH:" != *":${LIBTORCH_LIBDIR}:"* ]] && \
    #    export LD_LIBRARY_PATH="${LIBTORCH_LIBDIR}:${LD_LIBRARY_PATH}"
echo "Key env variables set"
echo "EDEPSIM_DIR=${EDEPSIM_DIR}"

## Custom Edits  
alias ls='ls --color=auto'
alias zlab='cd /cluster/tufts/wongjiradlabnu/zimani01'

export PS1='\[\e[1;32m\]Singularity~\[\e[1;34m\]/${PWD/*\//}\[\e[1;37m\]> \[\e[0m\]'