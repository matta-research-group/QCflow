import numpy as np
import rdkit
import itertools
import os
import os.path
import shutil
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors
import requests
from itertools import combinations
from typing import List


def write_gaussian(job_name, mol_name, smile, functional='B3LYP', basis_set='6-31G*', mol=None, 
                                                                                      torsion=None, 
                                                                                      conformer=None, 
                                                                                      old_chk=None):
    """
    job_name : The type of job run. Possible runs:
                Single Point neutral -> sp
                Optimisation neutral -> opt
                Torsional scan neutral → tor
                Optimisation neutral + Population analysis → pop_opt_n
                Single point anion → sp_a
                Single point cation → sp_c
                Optimisation anion → opt_a
                Optimisation cation → opt_c
                neutral charge, optimised anion geometry -> n_a_geo
                neutral charge, optimised cation geometry -> n_c_geo

    mol_name   : the name of the molecule

    smile      : The SMILE string of the molecule

    functional : Preset is B3LYP

    basis_set  : Preset is 6-31G*

    mol        : The rdkit molecule 
    
    torsion    : The getTorsion information in regards to the rotatable bond

    conformer  : The rdkit conformer of a molecule 

    old_chk     : The previous chk file used to guess geometry/MOs. Eg, 'mol_name_opt_321G.chk'
                
    """
    if (job_name=='sp'):

        old_chk = ' '
        # get atom names
        symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        # get x y z coords
        geometry = conformer
        torsion_data = f' \n'
        calculation = 'sp'
        
        if Descriptors.NumRadicalElectrons(mol) == 0:
            mult_chg = '0 1'
        if Descriptors.NumRadicalElectrons(mol) == 1:
            mult_chg = '0 2'
        if Descriptors.NumRadicalElectrons(mol) == 2:
            mult_chg = '0 3' 


    if (job_name=='opt'):

        if old_chk is not None:
            old_chk = old_chk # needed if starting from previous job
            
        else:
            old_chk = ' '
            symbols = [a.GetSymbol() for a in mol.GetAtoms()]
            # get x y z coords
            geometry = conformer
            
        torsion_data = f' \n'
        calculation = 'opt'
        if Descriptors.NumRadicalElectrons(mol) == 0:
            mult_chg = '0 1'
        if Descriptors.NumRadicalElectrons(mol) == 1:
            mult_chg = '0 2'
        if Descriptors.NumRadicalElectrons(mol) == 2:
            mult_chg = '0 3' 

    if (job_name=='tor'):

        old_chk = ' '
        # get atom names
        symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        # get x y z coords
        geometry=mol.GetConformers()[0]
        # unpack torsion terms and add 1 for fortran 1-based numbering
        a1, a2, a3, a4 = [t+1 for t in torsion]
        #default is to scan entire dihedral
        torsion_data = f'{a1} {a2} {a3} {a4} S 35 10.0\n'
        calculation = 'opt=modredundant'
        mult_chg = '0 1' # by default all molecules are neutral and singlets!

    if (job_name=='pop_opt_n'):

        if old_chk is not None:
            old_chk = old_chk # needed if starting from previous job
            calculation = 'Pop=Full, Geom=Checkpoint'
        else:
            old_chk = ' '
            # get atom names
            symbols = [a.GetSymbol() for a in mol.GetAtoms()]
            # get x y z coords
            geometry = conformer
            calculation = 'Pop=Full'
        torsion_data = f' \n'
        mult_chg = '0 1' # by default all molecules are neutral and singlets!

    if (job_name=='sp_a') or (job_name=='sp_c'):

        torsion_data = f' \n'
        # assumes old chk comes from opt neutral
        if old_chk is None:
            old_chk = f'%OldChk={mol_name}_opt.chk'
        calculation = 'SP, Geom=Checkpoint'
        #reads the geometry information from the checkpoint file

    if (job_name=='opt_c') or (job_name=='opt_a'):

        torsion_data = f' \n'
        if old_chk is None:
            old_chk = f'%OldChk={mol_name}_opt.chk'
        calculation = 'opt, Geom=Checkpoint'
        #reads the geometry information from the checkpoint file

    if (job_name=='sp_c') or (job_name=='opt_c'):
        mult_chg = '1 2' # This is a positive ion calc so multiplicity and charge change!
        
    if (job_name=='sp_a') or (job_name=='opt_a'):
        mult_chg = '-1 2' # This is a negative ion calc so multiplicity and charge change!
                                                                                          
    file_name = f'{mol_name}_{job_name}.com'
    #Includes information about basis set to allow for ramping
    chk_name = f'{mol_name}_{job_name}.chk'
    title = f'{mol_name} {job_name} Smile String: {smile}'

    with open(file_name, 'w') as file:
        file.write(f'{old_chk}\n')
        file.write(f'%Chk={chk_name}\n') #
        file.write(f'#p {functional}/{basis_set} {calculation}, nosymm\n') #
        file.write(' \n')#
        file.write(f'{title}\n')#
        file.write(' \n')#
        file.write(f'{mult_chg} \n')#
        # if no previous checkpoint file, get geometry
        if old_chk==' ':
            for atom,symbol in enumerate(symbols):
                p = geometry.GetAtomPosition(atom)
                # atom  x y z
                line = f' {symbol} {p.x:.5f} {p.y:.5f} {p.z:.5f} \n'
                file.write(line)
        file.write(' \n')
        file.write(torsion_data)
        file.write(' \n')
