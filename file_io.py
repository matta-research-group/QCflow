import MDAnalysis as mda
import itertools
import cclib
import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem
from itertools import combinations
from typing import List
from rdkit.Chem import AllChem, PeriodicTable


def submit_slurm_job(job_type, mol_name):
    """
    job_type : This is determined by what calculations you want to happen. The options are
                'torsional scan', 'neutral population', 'opt anion', 'ver anion', 'opt cation'
                and 'ver cation'

    mol_name : the name of the dimer from the dictionary e.g. if fragment 0 was attached to fragment 1
                    then the dimer name is 0_1

    """

    if (job_type=='torsional scan'):

        run_name = f'{mol_name}_torsion'

    if (job_type=='neutral population'):

        run_name = f'{mol_name}_neu_opt_pop'

    if (job_type=='opt anion'):

        run_name = f'{mol_name}_opt_anion'

    if (job_type=='ver anion'):

        run_name = f'{mol_name}_ver_anion'

    if (job_type=='opt cation'):

        run_name = f'{mol_name}_opt_cat'

    if (job_type=='ver cation'):

        run_name = f'{mol_name}_ver_cat'

    string = f'sbatch {run_name}.sh'

    process = subprocess.run(string,
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE, shell=True,check=True)
    return process.stdout
