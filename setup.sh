#!/bin/bash

if ! command -v conda &> /dev/null
then
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    	mkdir utils/miniconda3
    	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O utils/miniconda3/miniconda.sh
    elif [[ "$OSTYPE" == "darwin"* ]]; then
	    brew install libomp
	    mkdir utils/miniconda3
      curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o utils/miniconda3/miniconda.sh
    bash utils/miniconda3/miniconda.sh -b -u -p utils/miniconda3
    rm -rf utils/miniconda3/miniconda.sh
    CONDA_PATH="utils/miniconda3/bin/conda"
    CONDA_BASE="utils/miniconda3"
    CONDA_SH="utils/miniconda3/etc/profile.d/conda.sh"
    fi
else
    CONDA_PATH=$(which conda)
    CONDA_BASE=$(conda info --base)
    CONDA_SH=$CONDA_BASE/etc/profile.d/conda.sh
fi
$CONDA_PATH init bash
$CONDA_PATH config --set auto_activate_base false
echo -e "conda config --set channel_priority strict"
$CONDA_PATH config --set channel_priority strict
echo "conda env create -f scorch.yml python=3.6.9"
$CONDA_PATH env create -f scorch.yml python=3.6.9
