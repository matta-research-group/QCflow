![logo](https://github.com/matta-research-group/QCflow/blob/main/QCflow_logo_narrow.jpg?raw=true)

# QCflow

A cheminformatics -> quantum chemistry workflow toolkit leveraging [rdkit](https://github.com/rdkit/rdkit) and [cclib](https://github.com/cclib/cclib).

A typical workflow involves:

1. Generating a series of input molecules as `SMILES` - optionally, by combining different fragments into larger molecules/oligomers  

2. For each molecule:
    a. write `.com` Gaussian input files, `.sh` SLURM submission files 
    b. submit job (assuming you are working within a HPC)
    c. parse output file to submit further calculations or retrieve descriptors

3. Combine descriptors and results in a `pandas` dataframe format or similar for plotting / further analysis

## Supported QC codes 

Only Gaussian is supported at the moment, but we plan to add support for Psi4.


## Installation

```bash
git clone https://github.com/matta-research-group/QCflow.git
cd QCflow
# install requirements into new environment
conda env create -f QCflow.yml
conda activate QCflow
# install the QCflow package
pip install .
```

## Usage Examples

The `run_torsion` and `run_opt_neutral` functions contain example workflows that read in molecules from `.smi` files and submit quantum chemistry jobs for each molecule.

## Calculation settings

QCflow can prepare and submit input files for the following jobs: 
- Single point calculation -> `sp`
- Geometry optimisation -> `opt`
- Torsional scan neutral → `tor`
- Optimisation neutral/Population analysis → `pop_opt_n`
- Vertical anion → `ver_a`
- Vertical cation → `ver_c`
- Optimisation anion → `opt_a`
- Optimisation cation → `opt_c`
