#/usr/bin/bash

# detect the os and install dependencies
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        installcommand=$(command -v yum)
	if [ -z "$installcommand" ]; then
      		yes | sudo apt-get install libssl-dev
	else
      		yes | sudo yum install openssl libssl-dev
	fi
elif [[ "$OSTYPE" == "darwin"* ]]; then
        :
fi

# download python 3.6.9
echo "**************************************"
echo "Downloading Python 3.6.9..."
echo "**************************************"
curl -OL https://www.python.org/ftp/python/3.6.9/Python-3.6.9.tgz

# extract the downloaded python file
echo "**************************************"
echo "Extracting Python 3.6.9..."
echo "**************************************"
tar -xf Python-3.6.9.tgz && rm -f Python-3.6.9.tgz

# install python from source
echo "**************************************"
echo "Installing Python 3.6.9..."
echo "**************************************"
cd Python-3.6.9/ && ./configure -prefix=$PWD -enable-shared  && make && sudo make install && cd -

# make sure virtualenv is installed
echo "**************************************"
echo "Installing virtualenv..."
echo "**************************************"
python3 -m pip install virtualenv

# create scoring virtualenv
echo "**************************************"
echo "Installing virtual environment for scoring function..."
echo "**************************************"
virtualenv .scorch --python=Python-3.6.9/python

# activate the scoring virtualenv
echo "**************************************"
echo "Activating scoring function virtual environment"
echo "**************************************"
source .scorch/bin/activate

# install venv dependencies
echo "**************************************"
echo "Installing scoring function python dependencies"
echo "**************************************"
python -m pip install -r requirements.txt
