# SCORCH

[![License](https://img.shields.io/badge/License-MIT-blue)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Python-3.7.12-blue)](https://www.python.org/downloads/release/python-3712/)
![](https://img.shields.io/badge/OS-linux%20%7C%20OS%20X-blueviolet)
[![Preprocessing](https://img.shields.io/badge/preprocessing-MGLTools%201.5.6-brightgreen)](https://ccsb.scripps.edu/mgltools/1-5-6/)
[![Docking](https://img.shields.io/badge/docking-GWOVina%201.0-brightgreen)](https://doi.org/10.1111/cbdd.13764)

---

SCORCH is a fast scoring function based on a consensus of machine learning models. Scoring functions are used to evaluate poses of molecules obtained from molecular docking. SCORCH scores range from 0 to 1, with higher values indicating a higher probability of the molecule binding tightly to the receptor.

SCORCH uses `.pdbqt` files as input for the scoring, which is the format used by [Autodock](https://autodock.scripps.edu/), [Vina](https://dx.doi.org/10.1002/jcc.21334), and [GWOVina](https://cbbio.online/software/gwovina/index.html) docking software, among others. Additionally, this release contains an integrated pipeline to dock and score molecules in SMILES format using GWOVina.

SCORCH uses a variety of descriptors to characterize a docked pose, including [Binana 1.3](https://git.durrantlab.pitt.edu/jdurrant/binana/-/tree/1.3) and [ECIFs](https://github.com/DIFACQUIM/ECIF). The contributing models were trained on multiple docked poses for each ligand, labelled based on their RMSD to crystal structures. Training examples of non-binders were generated using [DeepCoy](https://github.com/fimrie/DeepCoy). SCORCH has used over 54,000 poses in its training. As a result, SCORCH avoids biases and provides improved accuracy to identify true binder molecules in virtual screening. Read more in our [publication]().


![](https://raw.githubusercontent.com/miles-mcgibbon/miles-mcgibbon/main/.github/images/pose_labels.gif)

---

# Installation

### **Currently SCORCH is only compatible with OS X and Linux; Docking and Scoring SMILES is only available on Linux systems**

Installation on linux and mac is achieved via [conda](https://docs.conda.io/en/latest/). The silent installation of miniconda (which will not affect any existing conda or python installations), SCORCH dependencies and SCORCH setup is all performed with the supplied setup bash script. This will additionally install MGLTools 1.5.6, and on Linux systems will also install GWOVina 1.0.

To install SCORCH:

```bash
# clone the GitHub repository
env GIT_SSL_NO_VERIFY=true git clone https://gitfront.io/r/mmcgibbon/cff5aca4804137fbf88d4c7357b32f0cd1c20ad1/SCORCH.git

# ensure setup.sh is executable
cd SCORCH
sudo chmod +x setup.sh

# execute the setup script
./setup.sh
```

# Receptor & Ligand Preparation

### **Please note MGLTools does not work with M1 Macs or Macs running Catalina OS**

The scoring function accepts `.pdbqt` receptor files and SMILES or `.pdbqt` ligand files as inputs. Any `.pdb` receptor files should be prepared with the supplied [MGLTools 1.5.6](https://ccsb.scripps.edu/mgltools/1-5-6/) using `prepare_receptor4.py` Python script as follows:

```bash

# preparing a receptor
./utils/MGLTools-1.5.6/bin/pythonsh \
./utils/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py \
-r path/to/receptor.pdb \
-A hydrogens \
-o path/to/save/receptor.pdbqt
-U nphs
```

SMILES ligands need no preprocessing.For pre-docked `.pdbqt` ligands, we recommend only scoring docking results in `.pdbqt` format (ideally from AutoDock or GWOVina). Scoring results from other docking software might be possible if converted to `.pdbqt`.

# Usage

To use the scoring function, the conda environment needs to be activated first:

```bash
conda activate scorch
```

The scoring function is supplied as the Python script `scorch.py`. Its main arguments are:

|Argument     |Value                                                                                     |Importance                  |
|-------------|------------------------------------------------------------------------------------------|----------------------------|
|`-receptor`    |Filepath to receptor file (pdbqt)                                                       |Essential                   |
|`-ligand `     |Filepath to ligand(s) (pdbqt or SMILES)                                                 |Essential                   |
|`-ref_lig`    |Filepath to example ligand in receptor binding site (mol, mol2, sdf, pdb or pdbqt)       |Essential for SMILES ligands (unless `-center` and `-range` supplied)|
|`-center`     | '[x, y, z]' coordinates of the center of the binding site for docking                       |Essential for SMILES ligands (unless `-ref_lig` supplied)|
|`-range`     | '[x, y, z]' axis lengths to define a box around `-center` coordinates for docking            |Essential for SMILES ligands (unless `-ref_lig` supplied)|
|`-out`         |Filepath for output csv (If not supplied, scores are written to stdout)                  |Optional (Default stdout)   |
|`-return_pose_scores` |If supplied, scoring values for individual poses in each ligand file are returned | Optional (Default False) |
|`-threads`     |Number of threads to use                                                                |Optional (Default 1)        |
|`-verbose`     |If supplied, progress bars and indicators are displayed while scoring                  |Optional (Default False)    |

For further details on arguments, run `python scorch.py -h`.

### Docking and Scoring SMILES Ligands Against a Receptor (**Linux Only**)

SCORCH includes a full pipeline to convert SMILES ligands to `.pdbqt` files using MGLTools 1.5.6, dock them using GWOVina, and score them with SCORCH:

```bash
python scorch.py \
-receptor /home/user/receptors/receptor.pdbqt \
-ligand /home/user/smiles/ligands.smi  \
-ref_lig /home/user/ligands/reference_ligand.pdb \
-out scoring_results.csv
```

The binding site for docking can be defined with a reference ligand as above with `-ref_lig`, or by supplying `-center` and `-range` values in the same way as for molecular docking software:

```bash
python scorch.py \
-receptor /home/user/receptors/receptor.pdbqt \
-ligand /home/user/smiles/ligands.smi  \
-center '[14, 19, -24]' \
-range '[14, 14, 14]' \
-out scoring_results.csv
```

The `-ligand` input should be supplied as `.smi`  or `.txt` file, with one SMILES ligand and an optional identifier per line, as in this example:

```bash
CCC(CC)O[C@@H]1C=C(C[C@H]([C@H]1NC(=O)C)O)C(=O)OCC 49817880
CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)[NH3+])C(=O)OCC 24848267
CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)N)C(=O)NOC 118722031
```

Docking settings can be changed by editing the `utils/params/dock_settings.json` file in the scoring function folder. See the function help for additional details.

### Scoring Already Docked Ligands Against a Receptor (**Linux & Mac**)

For scoring a single ligand - `/home/user/ligands/ligand.pdbqt` - against a single receptor - `/home/user/receptors/receptor.pdbqt`:

```bash
python scorch.py \
-receptor /home/user/receptors/receptor.pdbqt \
-ligand /home/user/ligands/ligand.pdbqt
```

For scoring all `.pdbqt` ligands in the directory - `/home/user/ligands/` - against a single receptor - `/home/user/receptors/receptor.pdbqt` - just supply the directory path to the `-ligand` argument:

```bash
python scorch.py \
-receptor /home/user/receptors/receptor.pdbqt \
-ligand /home/user/ligands/  \
-out scoring_results.csv
# this writes the output scores to a file
-out scoring_results.csv \
# this parallelises the scoring over 6 threads
-threads 6 \
# this displays scoring progress
-verbose
```

### Importing SCORCH as a Python module

The main function from `scorch.py` can be imported and used in other Python scripts. It takes a dictionary of parameters as inputs and returns a pandas dataframe of model scores identical to the normal scoring function output:

```python
from scorch import scoring, parse_module_args

input_parameters = {'ligand': '/path/to/ligand.pdbqt',
                  'receptor': '/path/to/receptor.pdbqt',
                  'threads': 4,
                  'verbose': False}

parsed_parameters = parse_module_args(input_parameters)

output = scoring(parsed_parameters)
```

# Output

Scores are output in `.csv` format. For example, scoring a single ligand `.pdbqt` file with 10 docked poses on a single receptor, using the `-return_pose_scores` flag, yielded the following output. Note that the output of SCORCH includes a measure of the prediction certainty.

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
