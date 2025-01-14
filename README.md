![logo](https://github.com/matta-research-group/QCflow/blob/main/QCflow_logo_narrow.jpg?raw=true)

# QCflow

[![QCflow stable](https://github.com/matta-research-group/QCflow/actions/workflows/run_test.yml/badge.svg?branch=qcflow-0.2)](https://github.com/matta-research-group/QCflow/actions/workflows/run_test.yml)
![Python](https://img.shields.io/badge/language-Python-blue.svg)
[![Documentation](https://img.shields.io/badge/Documentation-Online-brightgreen)](https://matta-research-group.github.io/QCflow/)
[![GitHub Last commit](https://img.shields.io/github/last-commit/matta-research-group/QCflow)](https://github.com/matta-research-group/QCflow/commits/qcflow-0.3)
[![GitHub stars](https://img.shields.io/github/stars/matta-research-group/QCflow)](https://github.com/matta-research-group/QCflow/stargazers)

A cheminformatics -> quantum chemistry workflow toolkit leveraging [rdkit](https://github.com/rdkit/rdkit) and [cclib](https://github.com/cclib/cclib).

A typical workflow involves:

1. Generating a series of input molecules as `SMILES` - optionally, by combining different fragments into larger molecules/oligomers  

2. For each molecule:
   
    a. write `.com` Gaussian input files, `.sh` SLURM submission files
    
    b. submit a job (assuming you are working within a HPC)
   
    c. parse output file to submit further calculations or retrieve descriptors

3. Combine descriptors and results in a [pandas](https://pandas.pydata.org/) dataframe format or similar for plotting / further analysis

## Supported QC codes 

Only [Gaussian](https://gaussian.com/man/) is supported at the moment, but we plan to add support for [Psi4](https://psicode.org/).

Psi4 is in the beta stages of support and can be found on the 'qcflow-psi4' branch.

When installing the branch with psi4 please use the following command onc you have installed QCflow:

```bash
conda update psi4
```

## Installation

```bash
git clone --single-branch --branch qcflow-psi4 https://github.com/matta-research-group/QCflow.git
cd QCflow
# install requirements into new environment
conda env create -f QCflow.yml
conda activate QCflow
# install the QCflow package
pip install .
```

## Usage Examples and Advice

The [run_calc](https://github.com/matta-research-group/QCflow/blob/qcflow-0.3/QCflow/run_gaussian.py#L50) function within `run_gaussian.py` is an example workflow that submit gaussian calculations to the [KCL CREATE HPC](https://www.kcl.ac.uk/research/facilities/hpc-digital-platforms). 

Users external to Kings College London will need to alter the slurm.py file to match their HPC submission requirements.

Altering the loaded modules found on line [55](https://github.com/matta-research-group/QCflow/blob/qcflow-0.3/QCflow/slurm.py#L55) within `slurm.py` file will allow you to submit jobs to your HPC.

```bash
file.write(f'module load gaussian_sse4/16-C-gcc-13.2.0 \n') ### KCL CREATE HPC

file.write(f'module load your_gaussian_module \n') ### Your HPC
```
Once the alteration has been made, just ```pip install .``` again and the package will update for your HPC.


## Calculation settings

QCflow can prepare and submit input files for the following jobs: 
- Single point calculation, neutral -> `sp`
- Single point calculation, anion → `sp_a`
- Single point calculation, cation → `sp_c`
- Single point calculation neutral charge, cationic geometry → `n_c_geo`
- Single point calculation neutral charge, anioinc geometry → `n_a_geo`
- Geometry optimisation, neutral -> `opt`
- Torsional scan, neutral → `tor`
- Optimisation anion → `opt_a`
- Optimisation cation → `opt_c`
- Optimisation neutral + Population analysis → `pop_opt_n`
- Single point Hirshfeld calculation → `sp_hirsh`

## Files

```bash
├── build
│   ├── bdist.macosx-10.13-x86_64
│   └── lib
│       └── QCflow
│           ├── energy_calculations.py
│           ├── find_torsion.py
│           ├── fragments.py
│           ├── __init__.py
│           ├── load_gaussian.py
│           ├── run_gaussian.py
│           ├── slurm.py
│           ├── testing_data.py
│           ├── torsion_parser.py
│           └── write_gaussian.py
├── LICENSE
├── QCflow
│   ├── energy_calculations.py
│   ├── find_torsion.py
│   ├── fragments.py
│   ├── future_functions
│   │   └── orbital_parse.py
│   ├── __init__.py
│   ├── load_gaussian.py
│   ├── __pycache__
│   │   ├── find_torsion.cpython-311.pyc
│   │   ├── find_torsion.cpython-39.pyc
│   │   ├── fragments.cpython-311.pyc
│   │   ├── fragments.cpython-38.pyc
│   │   ├── fragments.cpython-39.pyc
│   │   ├── __init__.cpython-311.pyc
│   │   ├── __init__.cpython-38.pyc
│   │   ├── __init__.cpython-39.pyc
│   │   ├── load_gaussian.cpython-311.pyc
│   │   ├── load_gaussian.cpython-39.pyc
│   │   ├── run_opt_neutral.cpython-311.pyc
│   │   ├── run_opt_neutral.cpython-39.pyc
│   │   ├── run_opt_set_dihedral.cpython-311.pyc
│   │   ├── slurm.cpython-311.pyc
│   │   ├── slurm.cpython-39.pyc
│   │   ├── torsion_parser.cpython-311.pyc
│   │   ├── torsion_parser.cpython-39.pyc
│   │   ├── torsion_run.cpython-311.pyc
│   │   ├── torsion_run.cpython-39.pyc
│   │   ├── write_gaussian.cpython-311.pyc
│   │   └── write_gaussian.cpython-39.pyc
│   ├── run_gaussian.py
│   ├── slurm.py
│   ├── testing_data.py
│   ├── torsion_parser.py
│   └── write_gaussian.py
├── qcflow.egg-info
│   ├── dependency_links.txt
│   ├── PKG-INFO
│   ├── SOURCES.txt
│   └── top_level.txt
├── QCflow_logo_narrow.jpg
├── QCflow.yml
├── README.md
├── setup.py
└── tests
    ├── b_18_v2
    │   ├── b_18_v2_n_a_geo.log
    │   ├── b_18_v2_n_c_geo.log
    │   ├── b_18_v2_opt_a.log
    │   ├── b_18_v2_opt_c.log
    │   ├── b_18_v2_opt.com
    │   ├── b_18_v2_opt.log
    │   ├── b_18_v2_opt.sh
    │   ├── b_18_v2_sp_a.log
    │   ├── b_18_v2_sp_c.log
    │   ├── b_18_v2_sp_hirsh.log
    │   └── b_18_v2_tor.log
    ├── example_dic.json
    ├── fake_dict.json
    ├── __pycache__
    │   ├── tests.cpython-312.pyc
    │   └── tests.cpython-312-pytest-7.4.4.pyc
    ├── test_dict.json
    ├── test.py
    └── torsion_1
        └── torsion_1_tor.log
```