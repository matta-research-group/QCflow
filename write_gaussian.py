import numpy as np
import rdkit
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


def write_gaussian(job_name, mol_name, smile, functional='wB97XD', basis_set='aug-cc-pVDZ', mol=None, torsion=None, conformer=None):
    """
    job_name : The type of job run. Possible runs:
                Torsional scan neutral → tor
                Optimisation neutral/Population analysis → pop_opt_n
                Vertical anion → ver_a
                Vertical cation → ver_c
                Optimisation anion → opt_a
                Optimisation cation → opt_c
                neutral optimised anion geometry -> n_a_geo
                neutral optimised cation geometry -> n_c_geo

    mol_name : the name of the oligomer as seen in the dictionary i.e. if melanin fragment (b) is combined
                with organic electronic fragment (1) in the v1 position will be b_1_v1

    smile : The SMILE string of the molecule

    functional : Preset is wB97XD

    basis_set : Preset is aug-cc-pVDZ

    mol : The rdkit string of the dimer

    torsion : The getTorsion information in regards to the rotatable bond

    conformer : The corrected rdkit coordinates of a dimer, only used after a torsional scan has been done and
                is the lowest energy coordinates
    """

    if (job_name=='tor'):

        old_chk = f' \n'
        # get atom names
        symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        # get x y z coords
        geometry=mol.GetConformers()[0]
        # unpack torsion terms and add 1 for fortran 1-based numbering
        a1, a2, a3, a4 = [t+1 for t in torsion]
        torsion_data = f'{a1} {a2} {a3} {a4} S 35 10.0\n'
        calculation = 'opt=modredundant'
        mult_chg = '0 1' # by default all molecules are neutral and singlets!

    if (job_name=='pop_opt_n'):

        # get atom names
        symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        # get x y z coords
        geometry = conformer
        old_chk = f' \n'
        torsion_data = f' \n'
        calculation = 'opt, Pop=Full'
        mult_chg = '0 1' # by default all molecules are neutral and singlets!

    if (job_name=='opt_a'):

        torsion_data = f' \n'
        old_chk = f'%OldChk={mol_name}_pop_opt_n.chk'
        calculation = 'opt, Guess=Read, Geom=AllCheckpoint'
        #reads the geometry information from the checkpoint file
        mult_chg = '-1 2' # This is a negative ion calc so multiplicty and charge change!

    if (job_name=='ver_a'):

        torsion_data = f' \n'
        old_chk = f'%OldChk={mol_name}_pop_opt_n.chk'
        calculation = 'SP, Guess=Read, Geom=AllCheckpoint'
        #reads the geometry information from the checkpoint file
        mult_chg = '-1 2' # This is a negative ion calc so multiplicty and charge change!

    if (job_name=='opt_c'):

        torsion_data = f' \n'
        old_chk = f'%OldChk={mol_name}_pop_opt_n.chk'
        calculation = 'opt, Guess=Read, Geom=AllCheckpoint'
        #reads the geometry information from the checkpoint file
        mult_chg = '1 2' # This is a positive ion so multiplicty and charge change!

    if (job_name=='ver_c'):

        torsion_data = f' \n'
        old_chk = f'%OldChk={mol_name}_pop_opt_n.chk'
        calculation = 'SP, Guess=Read, Geom=AllCheckpoint'
        #reads the geometry information from the checkpoint file
        mult_chg = '1 2' # This is a positive ion calc so multiplicty and charge change!

    if (job_name=='n_a_geo'):

        torsion_data = f' \n'
        old_chk = f'%OldChk={mol_name}_opt_a.chk'
        calculation = 'SP, Guess=Read, Geom=AllCheckpoint'
        #reads the geometry information from the checkpoint file
        mult_chg = '0 1'

    if (job_name=='n_c_geo'):
        torsion_data = f' \n'
        old_chk = f'%OldChk={mol_name}_opt_c.chk'
        calculation = 'SP, Guess=Read, Geom=AllCheckpoint'
        #reads the geometry information from the checkpoint file
        mult_chg = '0 1'



    file_name = f'{mol_name}_{job_name}.com'
    chk_name = f'{mol_name}_{job_name}_{basis_set}.chk'
    title = f'{mol_name} {job_name} Smile String: {smile}'

    with open(file_name, 'w') as file:
        file.write(f'{old_chk}\n')
        file.write(f'%nproc=4\n')
        file.write(f'%Chk={chk_name}\n') #
        file.write(f'#p {functional}/{basis_set} {calculation}, nosymm\n') #
        file.write(' \n')#
        file.write(f'{title}\n')#
        file.write(' \n')#
        file.write(f'{mult_chg} \n')#
        if (job_name=='tor') or (job_name=='pop_opt_n'):
            for atom,symbol in enumerate(symbols):
                p = geometry.GetAtomPosition(atom)
                # atom  x y z
                line = f' {symbol} {p.x:.5f} {p.y:.5f} {p.z:.5f} \n'
                file.write(line)
        file.write(' \n')
        file.write(torsion_data)
        file.write(' \n')
