#!/bin/bash
source ~/.bashrc
if ! command -v conda &> /dev/null
then
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    	sudo mkdir utils/miniconda3
    	sudo wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O utils/miniconda3/miniconda.sh
    elif [[ "$OSTYPE" == "darwin"* ]]; then
	    brew install libomp
	    sudo mkdir utils/miniconda3
      sudo curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o utils/miniconda3/miniconda.sh
    fi
    echo "bash utils etc..."
    sudo bash utils/miniconda3/miniconda.sh -b -u -p utils/miniconda3
    sudo rm -rf utils/miniconda3/miniconda.sh
    CONDA_PATH="utils/miniconda3/bin/conda"
    CONDA_BASE="utils/miniconda3"
    CONDA_SH="utils/miniconda3/etc/profile.d/conda.sh"
else
    CONDA_PATH=$(which conda)
    CONDA_BASE=$(conda info --base)
    CONDA_SH=$CONDA_BASE/etc/profile.d/conda.sh
fi
sudo $CONDA_PATH init bash
source $CONDA_SH
$CONDA_PATH config --set auto_activate_base false
conda config --set channel_priority strict
sudo conda env create -f scorch.yml python=3.6.9
