from fragments import *
from find_torsion import *
from write_gaussian import *
from slurm import *
from torsion_parser import *
from run_opt_neutral import *
import os
import json

def run_reorganisation_energy(mol_name, mol_smile, functional='wB97XD', basis_set='cc-pVDZ'):
    """
    Submits four calculations:  Vertical anion → ver_a
                                Vertical cation → ver_c
                                Optimisation anion → opt_a
                                Optimisation cation → opt_c

    mol_name : name of oligomer

    mol_smile : SMILE string of oligomer

    functional : preset is wB97XD

    basis_set : preset is cc-pVDZ

    Additional this function relies upon the checkpoint file of the _pop_opt_n
    and needs these checkpoint files need to be name <dimer_name>_pop_opt_n.chk
    """
    #goes into directory
    os.chdir(f'{mol_name}')

    for job_name in ['opt_c', 'opt_a', 'ver_c', 'ver_a']:
        #writes the gaussian file
        write_gaussian(job_name, mol_name, mol_smile, functional, basis_set)
        #writes the slurm script
        write_slurm(job_name, mol_name)
        #submits the job
        submit_slurm_job(job_name, mol_name)

    #goes back to previous directory
    os.chdir(os.path.dirname(os.getcwd()))
