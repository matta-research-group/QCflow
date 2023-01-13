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


def write_gaussian(job_type, mol_name, mol=None, torsion=None, conformer=None):
    """
    job_type : This is determined by what calculations you want to happen. The options are
                'torsional scan', 'neutral population', 'opt anion', 'ver anion', 'opt cation'
                and 'ver cation'

    mol_name : the name of the dimer from the dictionary e.g. if fragment 0 was attached to fragment 1
                    then the dimer name is 0_1

    mol : The rdkit string of the dimer

    torsion : The getTorsion information in regards to the rotatable bond

    conformer : The corrected rdkit coordinates of a dimer, only used after a torsional scan has been done and
                is the lowest energy coordinates
    """

    if (job_type=='torsional scan'):

        # get atom names
        symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        # get x y z coords
        geometry=mol.GetConformers()[0]
        # unpack torsion terms and add 1 for fortran 1-based numbering
        a1, a2, a3, a4 = [t+1 for t in torsion]
        torsion_data = f'{a1} {a2} {a3} {a4} S 35 10.0\n'
        file_name = f'{mol_name}_torsion.com'
        chk_name = f'{mol_name}_torsion.chk'
        title = f'scan dihedral {mol_name}'
        calculation = 'opt=modredundant'
        mult_chg = '0 1' # by default all dimers are neutral and singlets!

    if (job_type=='neutral opt'):

        # get atom names
        symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        # get x y z coords
        geometry = conformer
        torsion_data = f' \n'
        file_name = f'{mol_name}_neu_opt_pop.com'
        chk_name = f'{mol_name}_neu_opt_pop.chk'
        title = f'Neutral optimised with HOMO/LUMO {mol_name}'
        calculation = 'opt=modredundant, Pop=Full'
        mult_chg = '0 1' # by default all dimers are neutral and singlets!

    if (job_type=='opt anion'):

        torsion_data = f' \n'
        file_name = f'{mol_name}_opt_anion.com'
        chk_name = f'{mol_name}_neu_opt_pop.chk'
        title = f'Optimisation of anion dimer {mol_name}'
        calculation = 'opt=modredundant, Geom=Checkpoint'
        #reads the geomtry information from the checkpoint file
        mult_chg = '-1 2' # This is a negative dimer calc so multiplicty and charge change!

    if (job_type=='ver anion'):

        torsion_data = f' \n'
        file_name = f'{mol_name}_ver_anion.com'
        chk_name = f'{mol_name}_neu_opt_pop.chk'
        title = f'Veritcal of anion dimer {mol_name}'
        calculation = 'SP, Geom=Checkpoint'
        #reads the geomtry information from the checkpoint file
        mult_chg = '-1 2' # This is a negative dimer calc so multiplicty and charge change!

    if (job_type=='opt cation'):

        torsion_data = f' \n'
        file_name = f'{mol_name}_opt_cat.com'
        chk_name = f'{mol_name}_neu_opt_pop.chk'
        title = f'Optimisation of cation dimer {mol_name}'
        calculation = 'opt=modredundant, Geom=Checkpoint'
        #reads the geomtry information from the checkpoint file
        mult_chg = '1 2' # This is a negative dimer calc so multiplicty and charge change!

    if (job_type=='ver cation'):

        torsion_data = f' \n'
        file_name = f'{mol_name}_ver_cat.com'
        chk_name = f'{mol_name}_neu_opt_pop.chk'
        title = f'Vertical of cation dimer {mol_name}'
        calculation = 'SP, Geom=Checkpoint'
        #reads the geomtry information from the checkpoint file
        mult_chg = '1 2' # This is a negative dimer calc so multiplicty and charge change!


    with open(file_name, 'w') as file:
        file.write(f'%chk={chk_name}\n') #
        file.write(f'#p B3LYP/6-31G* {calculation}, nosymm\n') #
        file.write(' \n')#
        file.write(f'{title}\n')#
        file.write(' \n')#
        file.write(f'{mult_chg} \n')#
        if (job_type=='torsional scan') or (job_type=='neutral population'):
            for atom,symbol in enumerate(symbols):
                p = geometry.GetAtomPosition(atom)
                # atom  x y z
                line = f' {symbol} {p.x:.5f} {p.y:.5f} {p.z:.5f} \n'
                file.write(line)
        file.write(' \n')
        file.write(torsion_data)
        file.write(' \n')


def write_slurm(job_type, mol_name):

    '''
    Writes the slurm file for the relevent job type and with the mol name that
    is being submited

    job_type : This is determined by what calculations you want to happen. The options are
                'torsional scan', 'neutral population', 'opt anion', 'ver anion', 'opt cation'
                and 'ver cation'

    mol_name : the name of the dimer from the dictionary e.g. if fragment 0 was attached to fragment 1
                    then the dimer name is 0_1
    '''

    if (job_type=='torsional scan'):

        file_name = f'{mol_name}_torsion.sh'
        input_file = f'{mol_name}_torsion.com'

    if (job_type=='neutral population'):

        file_name = f'{mol_name}_neu_opt_pop.sh'
        input_file = f'{mol_name}_neu_opt_pop.com'

    if (job_type=='opt anion'):

        file_name = f'{mol_name}_opt_anion.sh'
        input_file = f'{mol_name}_opt_anion.com'

    if (job_type=='ver anion'):

        file_name = f'{mol_name}_ver_anion.sh'
        input_file = f'{mol_name}_ver_anion.com'

    if (job_type=='opt cation'):

        file_name = f'{mol_name}_opt_cat.sh'
        input_file = f'{mol_name}_opt_cat.com'

    if (job_type=='ver cation'):

        file_name = f'{mol_name}_ver_cat.sh'
        input_file = f'{mol_name}_ver_cat.com'


    title = f'#!/bin/bash --login'
    with open(file_name, 'w') as file:
        file.write(f'{title}\n')#
        file.write(f'#SBATCH -o g09_%J.out \n')
        file.write(f'#SBATCH -e g09_%J.err \n')#
        file.write(f'#SBATCH --job-name={mol_name} \n')
        file.write(f'#SBATCH -p cpu \n')
        file.write(f'#SBATCH --ntasks=10 \n')
        file.write(f'#SBATCH --nodes=1 \n')
        file.write(f'#SBATCH --tasks-per-node=10 \n')
        file.write(f'#SBATCH --mem-per-cpu=4000 \n')
        file.write(f'#SBATCH --time=48:00:00 \n')
        file.write(' \n')#
        file.write(f'INPUTFILE={input_file} \n')
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
