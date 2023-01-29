import MDAnalysis as mda
import itertools
import subprocess
import numpy as np
from rdkit import Chem
from itertools import combinations
from typing import List
from rdkit.Chem import AllChem, PeriodicTable

def write_slurm(job_name, mol_name):

    '''
    Writes the slurm file for the relevent job type and with the mol name that
    is being submited

    job_name : The type of job run. Possible runs
                Torsional scan neutral → tor
                Population analysis → pop_n
                Vertical anion → ver_a
                Vertical cation → ver_c
                Optimisation anion → opt_a
                Optimisation cation → opt_c
                Optimisation neutral → opt_n

    mol_name : the name of the dimer from the dictionary e.g. if fragment 0 was attached to fragment 1
                    then the dimer name is 0_1
    '''
    file_name = f'{mol_name}_{job_name}.sh'


    title = f'#!/bin/bash --login'
    with open(file_name, 'w') as file:
        file.write(f'{title}\n')#
        file.write(f'#SBATCH -o g09_%J.out \n')
        file.write(f'#SBATCH -e g09_%J.err \n')#
        file.write(f'#SBATCH --job-name={mol_name}_{job_name} \n')
        file.write(f'#SBATCH -p cpu \n')
        file.write(f'#SBATCH --ntasks=10 \n')
        file.write(f'#SBATCH --nodes=1 \n')
        file.write(f'#SBATCH --tasks-per-node=10 \n')
        file.write(f'#SBATCH --mem-per-cpu=4000 \n')
        file.write(f'#SBATCH --time=48:00:00 \n')
        file.write(' \n')#
        file.write(f'INPUTFILE={mol_name}_{job_name}.com \n')
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


def submit_slurm_job(job_name, mol_name):
    """
    job_name : The type of job run. Possible runs
                Torsional scan neutral → tor.com
                Population analysis → pop_n.com
                Vertical anion → ver_a.com
                Vertical cation → ver_c.com
                Optimisation anion → opt_a.com
                Optimisation cation → opt_c.com
                Optimisation neutral → opt_n.com

    mol_name : the name of the dimer from the dictionary e.g. if fragment 0 was attached to fragment 1
                    then the dimer name is 0_1

    """

    string = f'sbatch {mol_name}_{job_name}.sh'

    process = subprocess.run(string,
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE, shell=True,check=True)
    return process.stdout
