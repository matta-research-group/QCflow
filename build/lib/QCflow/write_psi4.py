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


def write_psi4(job_name, mol_name, smile, functional='b3lyp', basis_set='6-31g*', mol=None, conformer=None):
    """
    Generates a psi4 input file based on the provided parameters.
    
    Parameters
    ----------
    job_name (str): The type of job to run. Possible values:
        - 'sp': Single Point neutral
        - 'opt': Optimisation neutral                                                                                                                                        
    mol_name (str): The name of the molecule.
    smile (str): The SMILE string of the molecule.
    functional (str, optional): The functional to use. Default is 'B3LYP'.
    basis_set (str, optional): The basis set to use. Default is '6-31G*'.
    mol (rdkit.Chem.rdchem.Mol, optional): The RDKit embedded molecule object.
    conformer (rdkit.Chem.rdchem.Conformer, optional): The RDKit conformer of the molecule.
    
    Returns
    -------
         None: Writes the Gaussian input file to disk.

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

        old_chk = ' '
        symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        # get x y z coords
        geometry = conformer
        calculation = 'opt'
        torsion_data = f' \n'

        if Descriptors.NumRadicalElectrons(mol) == 0:
            mult_chg = '0 1'
        if Descriptors.NumRadicalElectrons(mol) == 1:
            mult_chg = '0 2'
        if Descriptors.NumRadicalElectrons(mol) == 2:
            mult_chg = '0 3' 
                                                                                          
    file_name = f'{mol_name}_{job_name}.py'
    #Includes information about basis set to allow for ramping
    title = f'{mol_name} {job_name} Smile String: {smile}'
    number_proc = 'set_num_threads(4)'
    memory_num = '4000 mb'

    with open(file_name, 'w') as file:
        file.write('import psi4\n')
        file.write(' \n')#
        file.write(f"psi4.set_memory('{memory_num}')\n")
        file.write('psi4.set_num_threads(4)\n')
        file.write(' \n')#
        file.write('mol_name = psi4.geometry("""\n')#
        file.write(f' {mult_chg} \n')#
        # if no previous checkpoint file, get geometry
        if old_chk==' ':
            for atom,symbol in enumerate(symbols):
                p = geometry.GetAtomPosition(atom)
                # atom  x y z
                line = f' {symbol} {p.x:.5f} {p.y:.5f} {p.z:.5f} \n'
                file.write(line)
        file.write('""") \n')
        file.write(' \n')
        file.write(f"psi4.set_options({{'basis': '{basis_set}', 'scf_type': 'df'}})\n")
        file.write(' \n')
        if (job_name=='sp'):
            file.write(f'energy({functional})')
        if (job_name=='opt'):
            file.write(f"energy, wfn = psi4.optimize('{functional}', return_wfn=True)  \n")
            file.write(' \n')
            file.write("optimized_geometry_xyz = wfn.molecule().to_string('xyz') \n")
            file.write(f"with open('{mol_name}_{job_name}.xyz', 'w') as xyz_file:\n")
            file.write("    xyz_file.write(optimized_geometry_xyz) \n")
            file.write(' \n')
            file.write(f"psi4.core.set_active_molecule(wfn.molecule()) \n")
            file.write(' \n')
            file.write(f"energy = psi4.energy('{functional}') \n")
            file.write('homo = wfn.epsilon_a_subset("AO", "ALL").np[-1] \n')
            file.write('lumo = wfn.epsilon_a_subset("AO", "ALL").np[0] \n')
            file.write('energy_gap = lumo - homo \n')
            file.write(' \n')
            file.write(f"with open('{mol_name}_{job_name}_energy_and_gap.txt', 'w') as file:\n")
            file.write('    file.write(f"Optimized energy: {energy:.6f} Hartree\\n") \n')
            file.write('    file.write(f"HOMO: {homo:.6f} Hartree\\n") \n')
            file.write('    file.write(f"LUMO: {lumo:.6f} Hartree\\n") \n')
            file.write('    file.write(f"Energy gap (HOMO-LUMO): {energy_gap:.6f} Hartree\n") \n')
        file.write(' \n')
