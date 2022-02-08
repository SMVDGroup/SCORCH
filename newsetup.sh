mkdir utils/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O utils/miniconda3/miniconda.sh
bash utils/miniconda3/miniconda.sh -b -u -p utils/miniconda3
rm -rf utils/miniconda3/miniconda.sh
utils/miniconda3/bin/conda init bash
utils/miniconda3/bin/conda init zsh
source ~/.bashrc
conda config --set auto_activate_base false
conda env create -f scorch.yml
conda activate scorch

