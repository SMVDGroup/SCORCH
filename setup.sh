#!/bin/bash

# define the directory SCORCH is cloned into
BASEDIR=$PWD

echo -e """
###############################################################
# installing development tools for compiling GWOVina & MGLTools
###############################################################
"""

# install dependencies for xgboost, GWOVina & MGLTools
if [[ "$OSTYPE" == "darwin"* ]]; then
    # dependencies for mac
    echo -e "\nDetected Max OS X!"
    brew install libomp boost make gcc

elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
    echo -e "\nDetected Linux OS!"
    # check if debian or redhat linux distro
    installcommand=$(command -v yum)

    # if yum doesn't exist, assume debian
  	if [ -z "$installcommand" ]; then
            echo -e "\nDetected Debian! Using apt-get to install packages..."
            yes | sudo apt-get update
            # install dependencies for debian with apt-get
        		yes | sudo apt-get install build-essential libboost-all-dev
  	else
            echo -e "\nDetected RedHat! Using yum to install packages..."
            # install dependencies for redhat with yum
            yes | sudo yum update
        		yes | sudo yum groupinstall 'Development Tools'
            yes | sudo yum install boost-devel
  	fi
fi

echo -e """
###############################################################
# verifying conda install or installing miniconda3 if not found
###############################################################
"""

# check if conda is installed, and install miniconda3 if not

# if conda is not a recognised command then download and install
if ! command -v conda &> /dev/null; then
    echo -e "\nNo conda found - installing..."
    # if linux then get linux version
    if [[ "$OSTYPE" == "linux"* ]]; then
    	sudo mkdir utils/miniconda3
    	sudo wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O utils/miniconda3/miniconda.sh
    # if mac then get mac version
    elif [[ "$OSTYPE" == "darwin"* ]]; then
	    sudo mkdir utils/miniconda3
      sudo curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o utils/miniconda3/miniconda.sh
    fi

    # install miniconda3
    sudo bash utils/miniconda3/miniconda.sh -b -u -p utils/miniconda3

    # remove the installer
    sudo rm -rf utils/miniconda3/miniconda.sh

    # define conda installation paths
    CONDA_PATH="utils/miniconda3/bin/conda"
    CONDA_BASE="utils/miniconda3"
    CONDA_SH="utils/miniconda3/etc/profile.d/conda.sh"
else
    echo -e "\nFound existing conda install!"
    # if conda not installed then find location of existing installation
    CONDA_PATH=$(which conda)
    CONDA_BASE=$(conda info --base)
    CONDA_SH=$CONDA_BASE/etc/profile.d/conda.sh
fi

echo -e """
###############################################################
# installing the SCORCH conda environment
###############################################################
"""
# source the bash files to enable conda command in the same session
source ~/.bashrc
source ~/.bash_profile

# initiate conda
$CONDA_PATH init bash

# source the conda shell script once initiated
source $CONDA_SH

# configure conda to install environment quickly and silently
$CONDA_PATH config --set auto_activate_base false
conda config --set channel_priority strict

# repair permissions from sudo installing conda
# (necessary as some users will get a permission error)
sudo chown -R $USER $CONDA_BASE
sudo chown -R $USER ~/.conda

# create the conda environment
conda env create -f scorch.yml python=3.6.9

echo -e """
###############################################################
# building GWOVina 1.0 in utils/
###############################################################
"""

# define which GWOVina release to install based on OS
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    PLATFORM="linux"
elif [[ "$OSTYPE" == "darwin"* ]]; then
    PLATFORM="mac"
fi

# build GWOVina 1.0 in the utils folder
cd utils && tar -xzvf gwovina-1.0.tar.gz
cd gwovina-1.0/build/$PLATFORM/release
sudo make -j2

# return to the base directory
cd $BASEDIR

echo -e """
###############################################################
# building MGLTools 1.5.6 in utils/
###############################################################
"""

# download the correct install tarball for user OS
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    PLATFORM='linux'
    MGLFOLDER="mgltools_x86_64Linux2_1.5.6"
    wget https://ccsb.scripps.edu/download/548/ -O utils/mgltools1-5-6.tar.gz --no-check-certificate

elif [[ "$OSTYPE" == "darwin"* ]]; then
    PLATFORM="mac"
    MGLFOLDER="mgltools_i86Darwin9_1.5.6"
    curl https://ccsb.scripps.edu/download/552/ -o utils/mgltools1-5-6.tar.gz -k
fi

# build MGLTools 1.5.6 in utils folder
cd utils && tar -xzvf mgltools1-5-6.tar.gz
mv $MGLFOLDER MGLTools-1.5.6
cd MGLTools-1.5.6 && sudo ./install.sh

# notify user setup is complete

echo -e "\nSCORCH setup complete!"
