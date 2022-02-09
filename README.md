# SCORCH

[![License](https://img.shields.io/badge/License-MIT-blue)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Python-3.6.9-blue)](https://www.python.org/downloads/release/python-369/)
![](https://img.shields.io/badge/OS-linux%20%7C%20OS%20X-blueviolet)
[![Preprocessing](https://img.shields.io/badge/preprocessing-MGLTools%201.5.6-brightgreen)](https://ccsb.scripps.edu/mgltools/1-5-6/)
[![Docking](https://img.shields.io/badge/docking-GWOVina%201.0-brightgreen)](https://doi.org/10.1111/cbdd.13764)

---

SCORCH is a fast scoring function based on a consensus of machine learning models. Scoring functions are used to evaluate poses of molecules obtained from molecular docking. SCORCH scores range from 0 to 1, with higher values indicating a higher probability of the molecule binding tightly to the receptor.

SCORCH uses `.pdbqt` files as input for the scoring, which is the format used by [Autodock](https://autodock.scripps.edu/), [Vina](https://dx.doi.org/10.1002/jcc.21334), and [GWOVina](https://cbbio.online/software/gwovina/index.html) docking software, among others. Additionally, this release contains an integrated pipeline to dock and score molecules in SMILES format using GWOVina.

SCORCH uses a variety of descriptors to characterize a docked pose, including [Binana 1.3](https://git.durrantlab.pitt.edu/jdurrant/binana/-/tree/1.3) and [ECIFs](https://github.com/DIFACQUIM/ECIF). The contributing models were trained on multiple docked poses for each ligand, labelled based on their RMSD to crystal structures. SCORCH has used over 54,000 poses in its training. As a result, SCORCH avoids biases and provides improved accuracy to identify true binder molecules in virtual screening. Read more in our [publication]().


![](https://raw.githubusercontent.com/miles-mcgibbon/miles-mcgibbon/main/.github/images/pose_labels.gif)

---

# Installation

Installation on linux and mac is achieved via [conda](https://docs.conda.io/en/latest/). The silent installation of miniconda (which will not affect any existing conda or python installations), SCORCH dependencies and SCORCH setup is all performed with the supplied setup bash script.

To install SCORCH:

```bash
# clone the GitHub repository
git clone https://github.com/miles-mcgibbon/SCORCH.git

# ensure setup.sh is executable
cd SCORCH
sudo chmod 755 setup.sh

# execute the setup script
sudo ./setup.sh
```

# Usage

## Function options

To use the scoring function, the conda environment needs to be activated first:

```bash
conda activate scorch
```

The scoring function is supplied as the Python script `scorch.py`. Its main arguments are:

|Argument       |Value                                                                                     |Importance                  |
|---------------|------------------------------------------------------------------------------------------|----------------------------|
|`-receptor`    |Filepath to receptor file (pdbqt)                                                       |Essential                   |
|`-ligand `     |Filepath to ligand(s) (pdbqt or SMILES)                                                 |Essential                   |
|`-center`      |List specifying the center point of the docking search box in the 3D frame of the receptor |Essential for SMILES ligands|
|`-box_size`    |List specifying the dimensions of the docking search box, in Angstrom                    |Essential for SMILES ligands|    
|`-out`         |Filepath for output csv (If not supplied, scores are written to stdout)                  |Optional (Default stdout)   |
|`-return_pose_scores` |If supplied, scoring values for individual poses in each ligand file are returned | Optional (Default False) |
|`-threads`     |Number of threads to use                                                                |Optional (Default 1)        |
|`-verbose`     |If supplied, progress bars and indicators are displayed while scoring                  |Optional (Default False)    |

Additional options are explained in the function help.


## Docking and Scoring Multiple SMILES Ligands Against a Single PDBQT Receptor

The module includes a full pipeline to convert SMILES ligands to `.pdbqt` files using MGLTools 1.5.6, dock them using GWOVina, and score them with SCORCH:

```bash
python scorch.py \
-receptor /home/user/receptors/receptor.pdbqt \
-ligand /home/user/smiles/ligands.smi  \
-center [4.53, 2.25, -7.28] \
-box_size [20., 20., 20.] \
-out scoring_results.csv
```

The `-ligand` input is provided as a `.smi` or `.txt` file, with one SMILES ligand and an optional identifier per line, as in this example:

```bash
CCC(CC)O[C@@H]1C=C(C[C@H]([C@H]1NC(=O)C)O)C(=O)OCC 49817880
CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)[NH3+])C(=O)OCC 24848267
CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)N)C(=O)NOC 118722031
```

The `-receptor` input is provided as a `.pdbqt` file. This can be prepared from a `.pdb` file with [MGLTools 1.5.6](https://ccsb.scripps.edu/mgltools/1-5-6/) using the `prepare_receptor4.py` script as follows:

```bash
# preparing a receptor
./utils/MGLTools-1.5.6/bin/pythonsh \
./utils/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py \
-l path/to/receptor.pdb \
-A hydrogens \
-o path/to/save/receptor.pdbqt
-U nphs
```

The docking search space can be specified using the `-box_size` and `-center` arguments. Alternatively, you can provide a reference ligand (in `.pdb`, `.pdbqt`, `.mol`, `.mol2` or `.sdf` format) on the target pocket of interest using the `-ref_lig` argument. The coordinates of this ligand will then be padded to define the search box.

Docking settings can be changed by editing the `utils/params/dock_settings.json` file in the scoring function folder. See the function help for additional details.


## Scoring already docked ligands

For scoring a single ligand - `/home/user/ligands/ligand.pdbqt` - on a single receptor - `/home/user/receptors/receptor.pdbqt`:

```bash
python scorch.py \
-receptor /home/user/receptors/receptor.pdbqt \
-ligand /home/user/ligands/ligand.pdbqt
```

For scoring all `.pdbqt` ligands in a directory - `/home/user/ligands/` - on a single receptor - `/home/user/receptors/receptor.pdbqt` - just supply the directory path to the `-ligand` argument:

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

The `-verbose` flag is used here with the `-out` flag, otherwise the progress indicators will be written to the results file.


## Importing SCORCH as a Python module

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

## Output

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
