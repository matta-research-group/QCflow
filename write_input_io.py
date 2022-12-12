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

def write_gaussian_scan_input(dimer, dimer_name, torsion):

    # get atom names
    symbols = [a.GetSymbol() for a in dimer.GetAtoms()]
    # get x y z coords
    geometry=dimer.GetConformers()[0]
    # unpack torsion terms and add 1 for fortran 1-based numbering
    a1, a2, a3, a4 = [t+1 for t in torsion]
    file_name = f'{dimer_name}.com'
    chk_name = f'{dimer_name}.chk'
    title = f'scan dihedral {dimer_name}'
    mult_chg = '0 1' # by default all dimers are neutral and singlets!
    with open(file_name, 'w') as file:
        file.write(f'%chk={chk_name}\n') #
        file.write('#p B3LYP/6-31G* opt=modredundant, nosymm\n') #
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
        file.write(f'{a1} {a2} {a3} {a4} S 35 10.0\n')
        file.write(' \n')

def write_slurm_input(dimer_name):


    file_name = f'{dimer_name}.sh'
    title = f'#!/bin/bash --login'
    with open(file_name, 'w') as file:
        file.write(f'{title}\n')#
        file.write(f'#SBATCH -o g09_%J.out \n')
        file.write(f'#SBATCH -e g09_%J.err \n')#
        file.write(f'#SBATCH --job-name={dimer_name}_job \n')
        file.write(f'#SBATCH -p cpu \n')
        file.write(f'#SBATCH --ntasks=10 \n')
        file.write(f'#SBATCH --nodes=1 \n')
        file.write(f'#SBATCH --tasks-per-node=10 \n')
        file.write(f'#SBATCH --mem-per-cpu=4000 \n')
        file.write(f'#SBATCH --time=48:00:00 \n')
        file.write(' \n')#
        file.write(f'INPUTFILE={dimer_name}.com \n')
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
