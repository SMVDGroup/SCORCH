#!/bin/bash

if [[ "$OSTYPE" == "darwin"* ]]; then
    brew install libomp
fi

if ! command -v conda &> /dev/null
then
    if [[ "$OSTYPE" == "linux"* ]]; then
    	sudo mkdir utils/miniconda3
    	sudo wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O utils/miniconda3/miniconda.sh
    elif [[ "$OSTYPE" == "darwin"* ]]; then
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
source ~/.bashrc
$CONDA_PATH init bash
source $CONDA_SH
$CONDA_PATH config --set auto_activate_base false
conda config --set channel_priority strict
sudo chown -R $USER $CONDA_BASE
sudo chown -R $USER ~/.conda
conda env create -f scorch.yml python=3.6.9

# setup GWOVina
echo "Unpacking GWOVina 1.0"
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    PLATFORM="linux"
elif [[ "$OSTYPE" == "darwin"* ]]; then
    PLATFORM="mac"
fi
cd utils && tar -xzvf gwovina-1.0.tar.gz
cd gwovina-1.0/build/$PLATFORM/release
echo "Building GWOVina 1.0"
sudo make -j2
cd --

echo -e "\nSCORCH setup complete!"

# build MGLTools
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    PLATFORM='linux'
    MGLFOLDER="mgltools_x86_64Linux2_1.5.6"
    wget https://ccsb.scripps.edu/download/548/ -O utils/mgltools1-5-6.tar.gz --no-check-certificate

elif [[ "$OSTYPE" == "darwin"* ]]; then
    PLATFORM="mac"
    MGLFOLDER="mgltools_i86Darwin9_1.5.6"
    curl https://ccsb.scripps.edu/download/552/ -o utils/mgltools1-5-6.tar.gz -k
fi
cd utils && tar -xzvf mgltools1-5-6.tar.gz

mv $MGLFOLDER MGLTools-1.5.6

cd MGLTools-1.5.6 && sudo ./install.sh
