# MLScore

[![License](https://img.shields.io/badge/License-MIT-blue)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Python-3.6.9-blue)](https://www.python.org/downloads/release/python-369/)
![](https://img.shields.io/badge/OS-linux%20%7C%20OS%20X-blueviolet)
[![Preprocessing](https://img.shields.io/badge/preprocessing-MGLTools%201.5.6-brightgreen)](https://ccsb.scripps.edu/mgltools/1-5-6/)
[![Docking](https://img.shields.io/badge/docking-GWOVina%201.0-brightgreen)](https://doi.org/10.1111/cbdd.13764)

---

MLScore is a scoring function based on a consensus of machine learning models. The contributing models were trained on multiple docked ligand poses labelled by RMSD, meaning only well docked strong binders are rewarded, and poorly docked strong binders and all weak binders are punished.

---

# Installation

Installation on linux and mac is achieved via [virtualenv](https://virtualenv.pypa.io/en/latest/). The installation of virtualenv, a local build of python 3.6.9, scoring function dependencies and scoring function setup is all performed with the supplied setup bash script.

To install MLScore:

```bash
# clone the GitHub repository
git clone https://github.com/miles-mcgibbon/MLSCORE.git

# ensure setup.sh is executable
cd MLSCORE
chmod 755 setup.sh

# execute the setup script
./setup.sh
```

# Usage

To use the scoring function, the virtual environment needs to be activated first:

```bash
source .scoring/bin/activate
```


## Overview

The scoring function itself is supplied as a python script `scoring.py`. Arguments are:

|Argument     |Value                                                                                     |Importance                  |
|-------------|------------------------------------------------------------------------------------------|----------------------------|
|`-receptor`    |Filepath to receptor file (pdbqt)                                                         |Essential                   |
|`-ligand `     |Filepath to ligand(s) (pdbqt or SMILES)                                                   |Essential                   |
|`-ref_lig`    |Filepath to example ligand in receptor binding site (pdb or pdbqt)                        |Essential for SMILES ligands|
|`-out`         |Filepath for output csv (If not supplied scores are written to stdout)                    |Optional (Default stdout)   |
|`-threads`     |Number of threads to use                                                                  |Optional (Default 1)        |
|`-num_networks`|Number of networks to average for mlpscore and widedeepscore model scores                 |Optional (Default 15)       |
|`-verbose`     |True if supplied, progress bars and indicators displayed while scoring                    |Optional (Default False)    |
|`-pose_1`      |True if supplied, only the first model in pdbqt ligand files is scored                    |Optional (Default False)    |
|`-detailed`    |True if supplied, outputs include scores from individual mlpscore and widedeepscore models|Optional (Default False)    |

## Examples
---
### Scoring a Single PDBQT Ligand Against a Single PDBQT Receptor
---

For scoring a single ligand - `/home/user/ligands/ligand.pdbqt` - against a single receptor - `/home/user/receptors/receptor.pdbqt`:

```bash
python scoring.py \
-receptor /home/user/receptors/receptor.pdbqt \
-ligand /home/user/ligands/ligand.pdbqt
```

This would print the results to stdout; they could instead be redirected to a file:

```bash
python scoring.py \
-receptor /home/user/receptors/receptor.pdbqt \
-ligand /home/user/ligands/ligand.pdbqt > scoring_results.csv
```

The `-out` argument could be used to achieve exactly the same result:

```bash
python scoring.py \
-receptor /home/user/receptors/receptor.pdbqt \
-ligand /home/user/ligands/ligand.pdbqt  \
-out scoring_results.csv
```

Scoring could also be sped up and monitored by using some optional flags:

```bash
python scoring.py \
-receptor /home/user/receptors/receptor.pdbqt \
-ligand /home/user/ligands/ligand.pdbqt  \
-out scoring_results.csv \
# this parallelises the scoring over 6 threads
-threads 6 \
# this displays scoring progress
-verbose
```

**Note - the `-verbose` flag must be used in conjunction with the `-out` flag, otherwise the progress indicators will be written to the results file**

### Scoring Multiple PDBQT Ligands Against a Single PDBQT Receptor
---

For scoring all ligands in the directory - `/home/user/ligands/` - against a single receptor - `/home/user/receptors/receptor.pdbqt` - just supply the directory filepath to the `-ligand` argument:


```bash
python scoring.py \
-receptor /home/user/receptors/receptor.pdbqt \
-ligand /home/user/ligands/  \
-out scoring_results.csv
```

### Scoring Multiple PDBQT Ligands Against Multiple PDBQT Receptors
---

The scoring function can also take a folder of folders as an argument to score multiple ligands against multiple receptors. For example, say you have a benchmarking dataset stored in the following way in the folder `dataset`:

```
dataset
└───TargetA
│   │   TargetA_receptor.pdbqt
│   │   TargetA_ligand.pdbqt
│
└───TargetB
    │   TargetB_receptor.pdbqt
    │   TargetB_ligand.pdbqt
```

As long as, in each subfolder, there is **one receptor file with the suffix receptor.pdbqt** and **one ligand file with the suffix ligand.pdbqt**, the scoring function can return scores for `TargetA_ligand.pdbqt` against `TargetA_receptor.pdbqt` and `TargetB_ligand.pdbqt` against `TargetB_receptor.pdbqt`. To do this, just supply the parent directory to both the `-ligand` and `-receptor` arguments:

```bash
python scoring.py \
-receptor /home/user/dataset/ \
-ligand /home/user/dataset/  \
-out scoring_results.csv
```

### Docking and Scoring Multiple SMILES Ligands Against a Single PDBQT Receptor
---

The scoring function also includes a full pipeline to convert SMILES ligands to 3D pdbqt files using [MGLTools 1.5.6](https://ccsb.scripps.edu/mgltools/1-5-6/), dock them using [GWOVina](https://doi.org/10.1111/cbdd.13764), and score them with MLScore.

Say you wish to determine whether some ligands could be potent influenza neuraminidase inhibitors. SMILES inputs should be supplied as .smi or .txt files, with one SMILES ligand and SMILE identifier per line separated by a space. For example, you have a .smi or .txt file of ligands as SMILES - `/home/user/smiles/influenza_inhibitors.smi` - which contains the following:

```bash
CCC(CC)O[C@@H]1C=C(C[C@H]([C@H]1NC(=O)C)O)C(=O)OCC 49817880
CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)[NH3+])C(=O)OCC 24848267
CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)N)C(=O)NOC 118722031
```

Alongside an .smi or .txt file, you need to supply a pdb or pdbqt ligand bound to the active site of the protein receptor file you are screening against as a reference for GWOVina docking using the `-ref_lig` argument:

```bash
python scoring.py \
-receptor /home/user/receptors/influenza_neuraminidase.pdbqt \
-ligand /home/user/smiles/influenza_inhibitors.smi  \
-ref_lig /home/user/ligands/oseltamivir.pdb
-out influenza_results.csv
```

Again, this can be done more quickly and informatively with additional user arguments:

```bash
python scoring.py \
-receptor /home/user/receptors/influenza_neuraminidase.pdbqt \
-ligand /home/user/smiles/influenza_inhibitors.smi  \
-ref_lig /home/user/ligands/oseltamivir.pdb
-out influenza_results.csv \
-threads 6 \
-verbose
```

By default, GWOVina docking settings are identical to those used in the publication. However, these can be changed by manually editing the `utils/params/dock_settings.json` file in the scoring function folder:

```json
{"gwovina_settings": {"exhaustiveness":3,
                     "num_wolves":1,
                     "num_modes":20,
                     "energy_range":4},
"padding":12}
```

The `padding` variable is a value in angstroms which is added to the center coordinate of the reference ligand to determine the docking site limits.

# Output Scores

Scores are output in csv format. For example, scoring a single ligand pdbqt containing 20 docked poses against a single receptor file would yield the following output:

|Receptor      |Ligand        |xgbscore_multi|mlpscore_multi_best_average|wdscore_multi_best_average|multi_consensus|multi_consensus_stdev|multi_consensus_range|
|--------------|--------------|--------------|---------------------------|--------------------------|---------------|---------------------|---------------------|
|receptor.pdbqt|ligand_pose_1 |0.9974552     |0.783138                   |0.7250332                 |0.83520883     |0.11715221           |0.27242196           |
|receptor.pdbqt|ligand_pose_2 |0.99668676    |0.77808344                 |0.73867756                |0.83781594     |0.113484696          |0.2580092            |
|receptor.pdbqt|ligand_pose_3 |0.99672925    |0.82014656                 |0.7103614                 |0.8424124      |0.11796457           |0.28636783           |
|receptor.pdbqt|ligand_pose_4 |0.9974542     |0.7051348                  |0.6221999                 |0.7749297      |0.16095018           |0.37525433           |
|receptor.pdbqt|ligand_pose_5 |0.92504144    |0.6327202                  |0.6123963                 |0.723386       |0.14283314           |0.31264514           |
|receptor.pdbqt|ligand_pose_6 |0.9980556     |0.8461044                  |0.74997765                |0.86471254     |0.10212856           |0.24807793           |
|receptor.pdbqt|ligand_pose_7 |0.9973322     |0.80237544                 |0.65094346                |0.81688374     |0.14178425           |0.34638876           |
|receptor.pdbqt|ligand_pose_8 |0.9990158     |0.8433915                  |0.757661                  |0.8666894      |0.099900395          |0.24135482           |
|receptor.pdbqt|ligand_pose_9 |0.48365158    |0.7629251                  |0.70412195                |0.65023285     |0.12021217           |0.2792735            |
|receptor.pdbqt|ligand_pose_10|0.004160531   |0.07061309                 |0.1636547                 |0.07947611     |0.06541413           |0.15949416           |
|receptor.pdbqt|ligand_pose_11|0.08189643    |0.6112932                  |0.47703212                |0.39007393     |0.22470208           |0.5293968            |
|receptor.pdbqt|ligand_pose_12|0.014417216   |0.2984496                  |0.3226113                 |0.21182604     |0.1399372            |0.3081941            |
|receptor.pdbqt|ligand_pose_13|0.016477646   |0.36441794                 |0.31394863                |0.23161475     |0.15351388           |0.3479403            |
|receptor.pdbqt|ligand_pose_14|0.010164384   |0.23590061                 |0.21350199                |0.15318899     |0.10154623           |0.22573623           |
|receptor.pdbqt|ligand_pose_15|0.01252496    |0.2042825                  |0.21431309                |0.14370686     |0.09284995           |0.20178813           |
|receptor.pdbqt|ligand_pose_16|0.9616012     |0.820614                   |0.71008056                |0.83076525     |0.10293345           |0.25152063           |
|receptor.pdbqt|ligand_pose_17|0.18216649    |0.26185253                 |0.18601812                |0.21001238     |0.03669023           |0.079686046          |
|receptor.pdbqt|ligand_pose_18|0.031303126   |0.51423347                 |0.4771736                 |0.3409034      |0.21944264           |0.48293033           |
|receptor.pdbqt|ligand_pose_19|0.0070598735  |0.22070605                 |0.2573203                 |0.1616954      |0.11036081           |0.25026044           |
|receptor.pdbqt|ligand_pose_20|0.4792845     |0.41491276                 |0.3693905                 |0.4211959      |0.045083493          |0.10989401           |
