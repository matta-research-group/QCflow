![logo](https://github.com/matta-research-group/QCflow/blob/main/QCflow_logo_narrow.jpg?raw=true)

# QCflow

[![QCflow stable](https://github.com/matta-research-group/QCflow/actions/workflows/run_test.yml/badge.svg?branch=qcflow-psi4)](https://github.com/matta-research-group/QCflow/actions/workflows/run_test.yml)
![Python](https://img.shields.io/badge/language-Python-blue.svg)
[![Documentation](https://img.shields.io/badge/Documentation-Online-brightgreen)](https://matta-research-group.github.io/QCflow/)
[![GitHub Last commit](https://img.shields.io/github/last-commit/matta-research-group/QCflow)](https://github.com/matta-research-group/QCflow/commits/qcflow-0.5/)
[![GitHub stars](https://img.shields.io/github/stars/matta-research-group/QCflow)](https://github.com/matta-research-group/QCflow/stargazers)

A cheminformatics -> quantum chemistry workflow toolkit leveraging [rdkit](https://github.com/rdkit/rdkit) and [cclib](https://github.com/cclib/cclib).

A typical workflow involves:

1. Generating a series of input molecules as `SMILES` - optionally, by combining different fragments into larger molecules/oligomers  

2. For each molecule a Gaussian 16 or Psi4 job is carried out:
   
    a. write `.com` Gaussian input files, `.sh` SLURM submission files

    b. submit a job (assuming you are working within a HPC)
   
    c. parse output file to submit further calculations or retrieve descriptors

    or

    a. write `.py` Psi4 input files, `.sh` SLURM submission files
    
    b. submit a job (assuming you are working within a HPC)
   
    c. parse output file to submit further calculations or retrieve descriptors

3. Combine descriptors and results in a [pandas](https://pandas.pydata.org/) dataframe format or similar for plotting / further analysis

## Installation

```bash
git clone https://github.com/matta-research-group/QCflow.git
cd QCflow
# install requirements into new environment
conda env create -f QCflow.yml
conda activate QCflow
# install the QCflow package
pip install .
conda update psi4
```

## Supported QC codes 

Only [Gaussian](https://gaussian.com/man/) and [Psi4](https://psicode.org/) are currently supported.

When installing QCflow check psi4 please by use the following command once you have installed QCflow to check psi4 is up to date:

```bash
conda update psi4
```

## Usage Examples and Advice

### Gaussian 16 Module

The [run_calc](https://github.com/matta-research-group/QCflow/blob/qcflow-0.3/QCflow/run_gaussian.py#L50) function within `run_gaussian.py` is an example workflow that submit gaussian calculations to the [KCL CREATE HPC](https://www.kcl.ac.uk/research/facilities/hpc-digital-platforms). 

Users external to King's College London will need to alter the slurm.py file to match their HPC submission requirements.

Altering the loaded modules found on line [59](https://github.com/matta-research-group/QCflow/blob/qcflow-0.3/QCflow/slurm.py#L59) within `slurm.py` file will allow you to submit jobs to your HPC.

```bash
file.write(f'module load gaussian_sse4/16-C-gcc-13.2.0 \n') ### KCL CREATE HPC

file.write(f'module load your_gaussian_module \n') ### Your HPC
```
Once the alteration has been made, just ```pip install .``` again and the package will update for your HPC.

### Psi4 Module

The [run_psi4](https://github.com/matta-research-group/QCflow/blob/qcflow-0.3/QCflow/run_psi4.py#L11) function within `run_psi4.py` is an example workflow that submit psi4 calculations to the [KCL CREATE HPC](https://www.kcl.ac.uk/research/facilities/hpc-digital-platforms).

Users external to King's College London will need to alter the slurm.py file to match their HPC submission requirements.

Altering the loaded modules found on line [117](https://github.com/matta-research-group/QCflow/blob/qcflow-0.3/QCflow/slurm.py#L117)

```bash
file.write(f'module load cuda/10.0.130-gcc-13.2.0 \n') ### KCL CREATE HPC

file.write(f'module load your_cuda_module \n') ### Your HPC
```

Once the alteration has been made, just ```pip install .``` again and the package will update for your HPC.

### SLURM Queue

The main queue within CREATE is `cpu` and this is denoted with the SLURM script `slurm.py` file. For external users to CREATE or users who wish to use a different queue, please change the `slurm.py` file to match your HPC submission requirements. The lines that need change are [48](https://github.com/matta-research-group/QCflow/blob/qcflow-0.3/QCflow/slurm.py#L48) and [109](https://github.com/matta-research-group/QCflow/blob/qcflow-0.3/QCflow/slurm.py#L109C9-L109C41) for Gaussian and Psi4 respectively.

```bash
file.write(f'#SBATCH -p cpu \n') #KCL CREATE HPC QUEUE

file.write(f'#SBATCH -p queue \n') #Your QUEUE
```

Once the alteration has been made, just ```pip install .``` again and the package will update for your HPC.

## Calculation settings - Gaussian 16

QCflow can prepare and submit input files for the following Gaussian 16 jobs: 
- Single point calculation, neutral → `sp`
- Single point calculation, anion → `sp_a`
- Single point calculation, cation → `sp_c`
- Single point calculation neutral charge, cationic geometry → `n_c_geo`
- Single point calculation neutral charge, anioinc geometry → `n_a_geo`
- Geometry optimisation, neutral → `opt`
- Torsional scan, neutral → `tor`
- Optimisation anion → `opt_a`
- Optimisation cation → `opt_c`
- Optimisation neutral + Population analysis → `pop_opt_n`
- Single point Hirshfeld calculation → `sp_hirsh`

## Calculation settings - Psi4

QCflow can prepare and submit input files for the following Psi4 jobs:
- Single point calculation, neutral → `sp`
- Geometry optimisation, neutral → `opt`
- Geometry optimisation, neutral from txt file geometry → `opt_pre_geom`
- Geometry optimisation cation (opt_c) and single of neutral charge, cation geometry (n_c_geo) → `cation`
- Geometry optimisation anion (opt_a) and single of neutral charge, anion geometry (n_a_geo) → `anion`
- Single point calculation, anion → `sp_a`
- Single point calculation, cation → `sp_c`

## Issues and Development

QCflow is under active development and we welcome any issues or suggestions. Please open an issue on the [GitHub repository](https://github.com/matta-research-group/QCflow/issues).

Additionally, we welcome any developments and improvements user may have and if you would like to contribute please open a pull request on the [GitHub repository](https://github.com/matta-research-group/QCflow/pulls).

## Files

```bash
.
├── build
│   ├── bdist.linux-x86_64
│   ├── bdist.macosx-11.0-arm64
│   └── lib
│       └── QCflow
│           ├── __init__.py
│           ├── energy_calculations.py
│           ├── find_torsion.py
│           ├── fragments.py
│           ├── load_gaussian.py
│           ├── run_gaussian.py
│           ├── run_psi4.py
│           ├── run_xTB.py
│           ├── sa_score.py
│           ├── slurm.py
│           ├── testing_data.py
│           ├── torsion_parser.py
│           ├── write_gaussian.py
│           ├── write_psi4.py
│           └── write_xTB.py
├── docs
│   ├── _build
│   │   ├── doctrees
│   │   │   ├── api.doctree
│   │   │   ├── environment.pickle
│   │   │   ├── getting_started.doctree
│   │   │   └── index.doctree
│   │   └── html
│   │       ├── _modules
│   │       │   ├── index.html
│   │       │   └── QCflow
│   │       │       ├── energy_calculations.html
│   │       │       ├── find_torsion.html
│   │       │       ├── fragments.html
│   │       │       ├── load_gaussian.html
│   │       │       ├── run_gaussian.html
│   │       │       ├── run_psi4.html
│   │       │       ├── sa_score.html
│   │       │       ├── slurm.html
│   │       │       ├── testing_data.html
│   │       │       ├── torsion_parser.html
│   │       │       ├── write_gaussian.html
│   │       │       └── write_psi4.html
│   │       ├── _sources
│   │       │   ├── api.rst.txt
│   │       │   ├── getting_started.rst.txt
│   │       │   └── index.rst.txt
│   │       ├── _static
│   │       │   ├── _sphinx_javascript_frameworks_compat.js
│   │       │   ├── basic.css
│   │       │   ├── css
│   │       │   │   ├── badge_only.css
│   │       │   │   ├── fonts
│   │       │   │   │   ├── fontawesome-webfont.eot
│   │       │   │   │   ├── fontawesome-webfont.svg
│   │       │   │   │   ├── fontawesome-webfont.ttf
│   │       │   │   │   ├── fontawesome-webfont.woff
│   │       │   │   │   ├── fontawesome-webfont.woff2
│   │       │   │   │   ├── lato-bold-italic.woff
│   │       │   │   │   ├── lato-bold-italic.woff2
│   │       │   │   │   ├── lato-bold.woff
│   │       │   │   │   ├── lato-bold.woff2
│   │       │   │   │   ├── lato-normal-italic.woff
│   │       │   │   │   ├── lato-normal-italic.woff2
│   │       │   │   │   ├── lato-normal.woff
│   │       │   │   │   ├── lato-normal.woff2
│   │       │   │   │   ├── Roboto-Slab-Bold.woff
│   │       │   │   │   ├── Roboto-Slab-Bold.woff2
│   │       │   │   │   ├── Roboto-Slab-Regular.woff
│   │       │   │   │   └── Roboto-Slab-Regular.woff2
│   │       │   │   └── theme.css
│   │       │   ├── doctools.js
│   │       │   ├── documentation_options.js
│   │       │   ├── file.png
│   │       │   ├── fonts
│   │       │   │   ├── Lato
│   │       │   │   │   ├── lato-bold.eot
│   │       │   │   │   ├── lato-bold.ttf
│   │       │   │   │   ├── lato-bold.woff
│   │       │   │   │   ├── lato-bold.woff2
│   │       │   │   │   ├── lato-bolditalic.eot
│   │       │   │   │   ├── lato-bolditalic.ttf
│   │       │   │   │   ├── lato-bolditalic.woff
│   │       │   │   │   ├── lato-bolditalic.woff2
│   │       │   │   │   ├── lato-italic.eot
│   │       │   │   │   ├── lato-italic.ttf
│   │       │   │   │   ├── lato-italic.woff
│   │       │   │   │   ├── lato-italic.woff2
│   │       │   │   │   ├── lato-regular.eot
│   │       │   │   │   ├── lato-regular.ttf
│   │       │   │   │   ├── lato-regular.woff
│   │       │   │   │   └── lato-regular.woff2
│   │       │   │   └── RobotoSlab
│   │       │   │       ├── roboto-slab-v7-bold.eot
│   │       │   │       ├── roboto-slab-v7-bold.ttf
│   │       │   │       ├── roboto-slab-v7-bold.woff
│   │       │   │       ├── roboto-slab-v7-bold.woff2
│   │       │   │       ├── roboto-slab-v7-regular.eot
│   │       │   │       ├── roboto-slab-v7-regular.ttf
│   │       │   │       ├── roboto-slab-v7-regular.woff
│   │       │   │       └── roboto-slab-v7-regular.woff2
│   │       │   ├── jquery.js
│   │       │   ├── js
│   │       │   │   ├── badge_only.js
│   │       │   │   ├── theme.js
│   │       │   │   └── versions.js
│   │       │   ├── language_data.js
│   │       │   ├── minus.png
│   │       │   ├── plus.png
│   │       │   ├── pygments.css
│   │       │   ├── README.md
│   │       │   ├── searchtools.js
│   │       │   └── sphinx_highlight.js
│   │       ├── api.html
│   │       ├── genindex.html
│   │       ├── getting_started.html
│   │       ├── index.html
│   │       ├── objects.inv
│   │       ├── py-modindex.html
│   │       ├── search.html
│   │       └── searchindex.js
│   ├── _static
│   │   └── README.md
│   ├── _templates
│   │   └── README.md
│   ├── api.rst
│   ├── conf.py
│   ├── getting_started.rst
│   ├── index.rst
│   ├── make.bat
│   ├── Makefile
│   ├── README.md
│   ├── requirements.yaml
│   └── timer.dat
├── Examples
│   ├── example_data_extraction_bulk.py
│   ├── example_data_extraction.ipynb
│   ├── example_g16_calculation_submission.py
│   ├── example_molecule_building.ipynb
│   ├── example_molecule_logs
│   │   ├── 18
│   │   │   └── 18_opt.log
│   │   ├── b
│   │   │   └── b_opt.log
│   │   └── b_18_single_v2
│   │       ├── b_18_single_v2_n_a_geo.log
│   │       ├── b_18_single_v2_n_c_geo.log
│   │       ├── b_18_single_v2_opt_a.log
│   │       ├── b_18_single_v2_opt_c.log
│   │       ├── b_18_single_v2_opt.log
│   │       ├── b_18_single_v2_sp_a.log
│   │       └── b_18_single_v2_sp_c.log
│   ├── example_psi4_opt_submission.py
│   ├── example_torsional_scan_script.py
│   └── FilteredSmilesFragments_props.csv
├── LICENSE
├── QCflow
│   ├── __init__.py
│   ├── __pycache__
│   │   ├── __init__.cpython-311.pyc
│   │   ├── __init__.cpython-313.pyc
│   │   ├── __init__.cpython-38.pyc
│   │   ├── __init__.cpython-39.pyc
│   │   ├── energy_calculations.cpython-313.pyc
│   │   ├── find_torsion.cpython-311.pyc
│   │   ├── find_torsion.cpython-313.pyc
│   │   ├── find_torsion.cpython-39.pyc
│   │   ├── fragments.cpython-311.pyc
│   │   ├── fragments.cpython-313.pyc
│   │   ├── fragments.cpython-38.pyc
│   │   ├── fragments.cpython-39.pyc
│   │   ├── load_gaussian.cpython-311.pyc
│   │   ├── load_gaussian.cpython-313.pyc
│   │   ├── load_gaussian.cpython-39.pyc
│   │   ├── run_gaussian.cpython-313.pyc
│   │   ├── run_opt_neutral.cpython-311.pyc
│   │   ├── run_opt_neutral.cpython-39.pyc
│   │   ├── run_opt_set_dihedral.cpython-311.pyc
│   │   ├── run_psi4.cpython-313.pyc
│   │   ├── run_xTB.cpython-313.pyc
│   │   ├── sa_score.cpython-313.pyc
│   │   ├── slurm.cpython-311.pyc
│   │   ├── slurm.cpython-313.pyc
│   │   ├── slurm.cpython-39.pyc
│   │   ├── testing_data.cpython-313.pyc
│   │   ├── torsion_parser.cpython-311.pyc
│   │   ├── torsion_parser.cpython-313.pyc
│   │   ├── torsion_parser.cpython-39.pyc
│   │   ├── torsion_run.cpython-311.pyc
│   │   ├── torsion_run.cpython-39.pyc
│   │   ├── write_gaussian.cpython-311.pyc
│   │   ├── write_gaussian.cpython-313.pyc
│   │   ├── write_gaussian.cpython-39.pyc
│   │   ├── write_psi4.cpython-313.pyc
│   │   └── write_xTB.cpython-313.pyc
│   ├── energy_calculations.py
│   ├── find_torsion.py
│   ├── fragments.py
│   ├── future_functions
│   │   └── orbital_parse.py
│   ├── load_gaussian.py
│   ├── run_gaussian.py
│   ├── run_psi4.py
│   ├── run_xTB.py
│   ├── sa_score.py
│   ├── slurm.py
│   ├── testing_data.py
│   ├── torsion_parser.py
│   ├── write_gaussian.py
│   ├── write_psi4.py
│   └── write_xTB.py
├── QCflow_logo_narrow.jpg
├── qcflow.egg-info
│   ├── dependency_links.txt
│   ├── PKG-INFO
│   ├── SOURCES.txt
│   └── top_level.txt
├── QCflow.yml
├── README.md
├── setup.py
└── tests
    ├── __pycache__
    │   ├── test.cpython-313-pytest-7.4.4.pyc
    │   ├── test.cpython-313-pytest-8.3.3.pyc
    │   ├── tests.cpython-312-pytest-7.4.4.pyc
    │   └── tests.cpython-312.pyc
    ├── 1
    │   ├── 1_cation.py
    │   ├── 1_opt.py
    │   └── 1_opt.sh
    ├── 12
    │   ├── 12_anion.err
    │   ├── 12_anion.out
    │   ├── 12_anion.py
    │   ├── 12_anion.sh
    │   ├── 12_cation.err
    │   ├── 12_cation.out
    │   ├── 12_cation.py
    │   ├── 12_cation.sh
    │   ├── 12_n_a_geo_energy_and_gap.txt
    │   ├── 12_n_c_geo_energy_and_gap.txt
    │   ├── 12_opt_a_energy_and_gap.txt
    │   ├── 12_opt_a.xyz
    │   ├── 12_opt_c_energy_and_gap.txt
    │   ├── 12_opt_c.xyz
    │   ├── 12_opt_energy_and_gap.txt
    │   ├── 12_opt.err
    │   ├── 12_opt.out
    │   ├── 12_opt.py
    │   ├── 12_opt.sh
    │   ├── 12_opt.xyz
    │   ├── 12_sp_a_energy_and_gap.txt
    │   ├── 12_sp_a.err
    │   ├── 12_sp_a.out
    │   ├── 12_sp_a.py
    │   ├── 12_sp_a.sh
    │   ├── 12_sp_c_energy_and_gap.txt
    │   ├── 12_sp_c.err
    │   ├── 12_sp_c.out
    │   ├── 12_sp_c.py
    │   ├── 12_sp_c.sh
    │   └── timer.dat
    ├── 145
    │   ├── 145_anion.err
    │   ├── 145_anion.out
    │   ├── 145_anion.py
    │   ├── 145_anion.sh
    │   ├── 145_cation.err
    │   ├── 145_cation.out
    │   ├── 145_cation.py
    │   ├── 145_cation.sh
    │   ├── 145_n_a_geo_energy_and_gap.txt
    │   ├── 145_n_c_geo_energy_and_gap.txt
    │   ├── 145_opt_a_energy_and_gap.txt
    │   ├── 145_opt_a.xyz
    │   ├── 145_opt_c_energy_and_gap.txt
    │   ├── 145_opt_c.xyz
    │   ├── 145_opt_energy_and_gap.txt
    │   ├── 145_opt.err
    │   ├── 145_opt.out
    │   ├── 145_opt.py
    │   ├── 145_opt.sh
    │   ├── 145_opt.xyz
    │   ├── 145_sp_a_energy_and_gap.txt
    │   ├── 145_sp_a.err
    │   ├── 145_sp_a.out
    │   ├── 145_sp_a.py
    │   ├── 145_sp_a.sh
    │   ├── 145_sp_c_energy_and_gap.txt
    │   ├── 145_sp_c.err
    │   ├── 145_sp_c.out
    │   ├── 145_sp_c.py
    │   ├── 145_sp_c.sh
    │   └── timer.dat
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
    ├── test_dict.json
    ├── test.py
    ├── timer.dat
    └── torsion_1
        └── torsion_1_tor.log
```