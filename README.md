# SCORCH

[![License](https://img.shields.io/badge/License-MIT-blue)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Python-3.7.12-blue)](https://www.python.org/downloads/release/python-3712/)
![](https://img.shields.io/badge/OS-linux%20%7C%20OS%20X-blueviolet)
[![Preprocessing](https://img.shields.io/badge/preprocessing-MGLTools%201.5.6-brightgreen)](https://ccsb.scripps.edu/mgltools/1-5-6/)
[![Docking](https://img.shields.io/badge/docking-GWOVina%201.0-brightgreen)](https://doi.org/10.1111/cbdd.13764)

---

SCORCH (Scoring COnsensus for RMSD-based Classification of Hits) is a fast scoring function based on a consensus of machine learning models. Scoring functions are used to evaluate poses of molecules obtained from molecular docking. SCORCH scores range from 0 to 1, with higher values indicating a higher probability of the molecule binding tightly to the receptor.

SCORCH uses `.pdbqt` files as input for the scoring, which is the format used by [Autodock](https://autodock.scripps.edu/), [Vina](https://dx.doi.org/10.1002/jcc.21334), and [GWOVina](https://cbbio.online/software/gwovina/index.html) docking software, among others. Additionally, this release contains an integrated pipeline to dock and score molecules in SMILES format using GWOVina.

SCORCH uses a variety of descriptors to characterize a docked pose, including [Binana 1.3](https://git.durrantlab.pitt.edu/jdurrant/binana/-/tree/1.3) and [ECIFs](https://github.com/DIFACQUIM/ECIF). The contributing models were trained on multiple docked poses for each ligand, labelled based on their RMSD to crystal structures. Training examples of non-binders were generated using [DeepCoy](https://github.com/fimrie/DeepCoy). SCORCH has used over 54,000 poses in its training. As a result, SCORCH avoids biases and provides improved accuracy to identify true binder molecules in virtual screening. Read more in our [publication](https://doi.org/10.1016/j.jare.2022.07.001).


![](https://raw.githubusercontent.com/miles-mcgibbon/miles-mcgibbon/main/.github/images/pose_labels.gif)

---

# Installation

## Linux

Installation on Linux is achieved via [conda](https://docs.conda.io/en/latest/). The supplied setup bash script installs SCORCH and dependencies in a conda environment, MGLTools 1.5.6 and GWOVina 1.0. If you do not have conda installed, the setup bash script installs it for you silently in the SCORCH directory.

To install SCORCH on Linux:

```bash
# clone the GitHub repository
git clone https://github.com/SMVDGroup/SCORCH.git

# ensure setup.sh is executable
cd SCORCH
sudo chmod +x setup.sh

# execute the setup script
./setup.sh
```

## Mac OS

Due to MGLTools-1.5.6 conflicts with newer Mac OS versions and M1 chips, SCORCH and SCORCH's SMILES docking and scoring pipeline can only be run on Mac OS systems via [Docker](https://www.docker.com/) which can be installed [here](https://docs.docker.com/get-docker/). Docker can also be used to run SCORCH on Linux systems in case of problems with the installation script.

To set up SCORCH once Docker is installed:

```bash
# download the scorch docker image
docker pull scorchml/scorch:v1.0

# run the scorch docker image
docker run -i -t scorchml/scorch:v1.0 /bin/bash

# once inside the docker image
cd home/SCORCH
```

# Receptor & Ligand Preparation

The scoring function accepts `.pdbqt` receptor files and SMILES or `.pdbqt` ligand files as inputs. Any `.pdb` receptor files should be prepared with the supplied [MGLTools 1.5.6](https://ccsb.scripps.edu/mgltools/1-5-6/) using `prepare_receptor4.py` Python script as follows:

```bash

# preparing a receptor
./utils/MGLTools-1.5.6/bin/pythonsh \
./utils/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py \
-r examples/predocked_1a0q/1a0q_receptor.pdb \
-A hydrogens \
-o examples/predocked_1a0q/1a0q_receptor_out.pdbqt \
-U nphs
```

SMILES ligands need no preprocessing. For pre-docked `.pdbqt` ligands, we recommend only scoring docking results in `.pdbqt` format (ideally from AutoDock or GWOVina). Scoring results from other docking software might be possible if converted to `.pdbqt`.

# Usage

**To use the scoring function on both Linux and Mac OS, the conda environment needs to be activated first:**

```bash
conda activate scorch
```

The scoring function is supplied as the Python script `scorch.py`. Its main arguments are:

|Argument     |Value                                                                                     |Importance                  |
|-------------|------------------------------------------------------------------------------------------|----------------------------|
|`-r, --receptor`    |Filepath to receptor file (pdbqt)                                                       |Essential                   |
|`-l, --ligand `     |Filepath to ligand(s) (pdbqt or SMILES)                                                 |Essential                   |
|`-rl, --ref_lig`    |Filepath to example ligand in receptor binding site (mol, mol2, sdf, pdb or pdbqt)       |Essential for SMILES ligands (unless `--center` and `--range` supplied)|
|`-c, --center`     | '[x, y, z]' coordinates of the center of the binding site for docking                       |Essential for SMILES ligands (unless `--ref_lig` supplied)|
|`-ra, --range`     | '[x, y, z]' axis lengths to define a box around `--center` coordinates for docking            |Essential for SMILES ligands (unless `--ref_lig` supplied)|
|`-o, --out`         |Filepath for output csv (If not supplied, scores are written to stdout)                  |Optional (Default stdout)   |
|`-p, --return_pose_scores` |If supplied, scoring values for individual poses in each ligand file are returned | Optional (Default False) |
|`-t, --threads`     |Number of threads to use                                                                |Optional (Default 1)        |
|`-v, --verbose`     |If supplied, progress bars and indicators are displayed while scoring                  |Optional (Default False)    |

For further details on arguments, run `python scorch.py --h`.

### Docking and Scoring SMILES Ligands Against a Receptor

SCORCH includes a full pipeline to convert SMILES ligands to `.pdbqt` files using MGLTools 1.5.6, dock them using GWOVina, and score them with SCORCH:

```bash
python scorch.py \
--receptor examples/smiles_REPTIN/reptin_receptor.pdbqt \
--ligand examples/smiles_REPTIN/reptin_smiles.smi  \
--ref_lig examples/smiles_REPTIN/reptin_ref_lig.pdbqt \
--out scoring_results.csv
```

The binding site for docking can be defined with a reference ligand as above with `--ref_lig`, or by supplying `--center` and `--range` values in the same way as for molecular docking software:

```bash
python scorch.py \
--receptor examples/smiles_REPTIN/reptin_receptor.pdbqt \
--ligand examples/smiles_REPTIN/reptin_smiles.smi  \
--center '[23.981,-42.667,67.156]' \
--range '[16.5,14.25,16.5]' \
--out scoring_results.csv
```

The `--ligand` input should be supplied as `.smi`  or `.txt` file, with one SMILES ligand and an optional identifier per line, as in this example:

```bash
CCC(CC)O[C@@H]1C=C(C[C@H]([C@H]1NC(=O)C)O)C(=O)OCC 49817880
CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)[NH3+])C(=O)OCC 24848267
CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)N)C(=O)NOC 118722031
```

Docking settings can be changed by editing the `utils/params/dock_settings.json` file in the scoring function folder. See the function help for additional details.

### Scoring Already Docked Ligands Against a Receptor

For scoring a single ligand - `examples/predocked_1a0q/ligands/1a0q_docked_ligand.pdbqt` - against a single receptor - `examples/predocked_1a0q/1a0q_receptor.pdbqt`:

```bash
python scorch.py \
--receptor examples/predocked_1a0q/1a0q_receptor.pdbqt \
--ligand examples/predocked_1a0q/ligands/1a0q_docked_ligand.pdbqt
```

For scoring all `.pdbqt` ligands in the directory - `examples/predocked_1a0q/ligands/` - against a single receptor - `examples/predocked_1a0q/1a0q_receptor.pdbqt` - just supply the directory path to the `--ligand` argument:

```bash
python scorch.py \
--receptor examples/predocked_1a0q/1a0q_receptor.pdbqt \
--ligand examples/predocked_1a0q/ligands/ \
# this writes the output scores to a file
--out scoring_results.csv \
# this parallelises the scoring over 6 threads
--threads 6 \
# this displays scoring progress
--verbose
```

### Importing SCORCH as a Python module

The main function from `scorch.py` can be imported and used in other Python scripts. It takes a dictionary of parameters as inputs and returns a pandas dataframe of model scores identical to the normal scoring function output:

```python
from scorch import scoring, parse_module_args

input_parameters = {'ligand': 'examples/predocked_1a0q/ligands/',
                  'receptor': 'examples/predocked_1a0q/1a0q_receptor.pdbqt',
                  'threads': 4,
                  'verbose': False}

parsed_parameters = parse_module_args(input_parameters)

output = scoring(parsed_parameters)

print(output)
```

# Output

Scores are output in `.csv` format. For example, scoring a single ligand `.pdbqt` file with 10 docked poses on a single receptor, using the `--return_pose_scores` flag, yielded the following output. Note that the output of SCORCH includes a measure of the prediction certainty.

|Receptor      |Ligand        |SCORCH_pose_score|SCORCH_certainty|Ligand_ID|SCORCH_score|best_pose|
|--------------|--------------|-----------------|----------------|---------|------------|---------|
|receptor.pdbqt|ligand_pose_1 |0.83521          |0.75148         |ligand   |0.86669     |0        |
|receptor.pdbqt|ligand_pose_2 |0.83782          |0.75926         |ligand   |0.86669     |0        |
|receptor.pdbqt|ligand_pose_3 |0.84241          |0.74976         |ligand   |0.86669     |0        |
|receptor.pdbqt|ligand_pose_4 |0.77493          |0.65857         |ligand   |0.86669     |0        |
|receptor.pdbqt|ligand_pose_5 |0.72339          |0.697           |ligand   |0.86669     |0        |
|receptor.pdbqt|ligand_pose_6 |0.86471          |0.78335         |ligand   |0.86669     |0        |
|receptor.pdbqt|ligand_pose_7 |0.81688          |0.69923         |ligand   |0.86669     |0        |
|receptor.pdbqt|ligand_pose_8 |0.86669          |0.78808         |ligand   |0.86669     |1        |
|receptor.pdbqt|ligand_pose_9 |0.65023          |0.74499         |ligand   |0.86669     |0        |
|receptor.pdbqt|ligand_pose_10|0.07948          |0.86123         |ligand   |0.86669     |0        |
