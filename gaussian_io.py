import numpy as np
import matplotlib.pyplot as plt
import rdkit
import cclib
import itertools
import os
import os.path
import shutil
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import requests
from itertools import combinations
from typing import List

def write_gaussian(job_type, mol_name, dimer=None, torsion=None, conformer=None):
    """
    Writes guassian files for specified job

    job_type : This is determined by what calculations you want to happen. The options are
                'torsional scan', 'neutral population', 'opt anion', 'ver anion', 'opt cation'
                and 'ver cation'

    mol_name : the name of the molecule from the dictionary e.g. if fragment 0 was attached to fragment 1
                    then the dimer name is 0-1

    mol : The rdkit Molecule 

    torsion : The getTorsion information in regards to the rotatable bond

    conformer : contains the rdkit Molecule with its update geometry
    """

    if (job_type=='torsional scan'):

        # get atom names
        symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        # get x y z coords
        geometry=dimer.GetConformers()[0]
        # unpack torsion terms and add 1 for fortran 1-based numbering
        a1, a2, a3, a4 = [t+1 for t in torsion]
        torsion_data = f'{a1} {a2} {a3} {a4} S 35 10.0\n'
        file_name = f'{mol_name}_torsion.com'
        chk_name = f'{mol_name}_torsion.chk'
        title = f'scan dihedral {mol_name}'
        calculation = 'opt=modredundant' # the only time we need modredundant
        mult_chg = '0 1' # by default all dimers are neutral and singlets!

    if (job_type=='neutral opt'):

        # get atom names
        symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        # get x y z coords
        geometry = conformer
        torsion_data = f' \n'
        file_name = f'{dimer_name}_neu_opt_pop.com'
        chk_name = f'{dimer_name}_neu_opt_pop.chk'
        title = f'Neutral optimised with HOMO/LUMO {dimer_name}'
        calculation = 'opt, Pop=Full'
        mult_chg = '0 1' # by default all dimers are neutral and singlets!

    if (job_type=='opt anion'):

        # get atom names
        symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        geometry = conformer
        torsion_data = f' \n'
        file_name = f'{mol_name}_opt_anion.com'
        chk_name = f'{mol_name}_opt_anion.chk'
        title = f'Optimisation of anion dimer {mol_name}'
        mult_chg = '-1 2' # This is a negative dimer calc so multiplicity and charge change!

    if (job_type=='ver anion'):

        # get atom names
        symbols = [a.GetSymbol() for a in dimer.GetAtoms()]
        geometry = comformer
        torsion_data = f' \n'
        file_name = f'{mol_name}_ver_anion.com'
        chk_name = f'{mol_name}_ver_anion.chk'
        title = f'Veritcal of anion dimer {mol_name}'
        calculation = 'sp'
        mult_chg = '-1 2' # This is a negative dimer calc so multiplicty and charge change!

    if (job_type=='opt cation'):

        # get atom names
        symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        geometry = comformer
        torsion_data = f' \n'
        file_name = f'{mol_name}_opt_cat.com'
        chk_name = f'{mol_name}_opt_cat.chk'
        title = f'Optimisation of cation dimer {mol_name}'
        calculation = 'opt'
        mult_chg = '1 2' # This is a negative dimer calc so multiplicty and charge change!

    if (job_type=='ver cation'):

        # get atom names
        symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        geometry = conformer
        torsion_data = f' \n'
        file_name = f'{mol_name}_ver_cat.com'
        chk_name = f'{mol_name}_ver_cat.chk'
        title = f'Vertical of cation dimer {mol_name}'
        calculation = 'sp'
        mult_chg = '1 2' # This is a negative dimer calc so multiplicty and charge change!


    with open(file_name, 'w') as file:
        file.write(f'%chk={chk_name}\n') #
        file.write(f'#p B3LYP/6-31G* {calculation}, nosymm\n') #
        file.write(' \n')#
        file.write(f'{title}\n')#
        file.write(' \n')#
        file.write(f'{mult_chg} \n')#
        for atom,symbol in enumerate(symbols):
            p = geometry.GetAtomPosition(atom)
            # atom  x y z
            line = f' {symbol} {p.x:.5f} {p.y:.5f} {p.z:.5f} \n'
            file.write(line)
        file.write(' \n')
        file.write(torsion_data)
        file.write(' \n')

def write_slurm(job_name):
    '''
    Provide the function with the name of the job and it produces a SLURM file for that job

    job_name : e.g. 0_1_ver_anion
    '''


    file_name = f'{job_name}.sh'
    title = f'#!/bin/bash --login'
    with open(file_name, 'w') as file:
        file.write(f'{title}\n')#
        file.write(f'#SBATCH -o g09_%J.out \n')
        file.write(f'#SBATCH -e g09_%J.err \n')#
        file.write(f'#SBATCH --job-name={job_name} \n')
        file.write(f'#SBATCH -p cpu \n')
        file.write(f'#SBATCH --ntasks=10 \n')
        file.write(f'#SBATCH --nodes=1 \n')
        file.write(f'#SBATCH --tasks-per-node=10 \n')
        file.write(f'#SBATCH --mem-per-cpu=4000 \n')
        file.write(f'#SBATCH --time=48:00:00 \n')
        file.write(' \n')#
        file.write(f'INPUTFILE={job_name}.com \n')
        file.write(f'OUTPUTFILE="\$(basename "\$INPUTFILE" .com).log" \n')
        file.write(' \n')
        file.write(f'module purge \n')
        file.write(f'module load test_switch_kcl/1.0.0-gcc-9.4.0 \n')
        file.write(f'module load gaussian_sse4_kcl/09-E-gcc-9.4.0 \n')
        file.write(f'export GOMP_CPU_AFFINITY=$SGE_BINDING \n')
        file.write(f'export KMP_AFFINITY="explicit,proclist=$SGE_BINDING,verbose" \n')
        file.write(f'#source $g09root/bsd/g09.login \n')
        file.write(' \n')
        file.write(f'echo "G09 job \$SLURM_JOBID" \n')
        file.write(f'echo "INPUT \$INPUTFILE" \n')
        file.write(f'echo "OUTPUT \$OUTPUTFILE" \n')
        file.write(f'echo "Running \$SLURM_NTASKS on \$SLURM_JOB_NODELIST" \n')
        file.write(' \n')
        file.write(f'#Execution Line \n')
        file.write(f'g09 "$INPUTFILE" > "$OUTPUTFILE" \n')
