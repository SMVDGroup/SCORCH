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
    utils/miniconda3/bin/conda init bash
    utils/miniconda3/bin/conda init zsh
    source ~/.bashrc
    conda deactivate
    conda config --set auto_activate_base false
    fi
fi
conda config --set channel_priority strict
conda env create -f scorch.yml python=3.6.9
conda activate scorch
