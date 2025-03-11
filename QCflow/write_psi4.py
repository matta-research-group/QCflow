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
    functional (str, optional): The functional to use. Default is 'b3lyp'.
    basis_set (str, optional): The basis set to use. Default is '6-31g'.
    mol (rdkit.Chem.rdchem.Mol, optional): The RDKit embedded molecule object.
    conformer (rdkit.Chem.rdchem.Conformer, optional): The RDKit conformer of the molecule.
    
    Returns
    -------
         None: Writes the psi4 input file to disk.

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
            file.write(f"energy, wfn = psi4.energy('{functional}', return_wfn=True)  \n")
            file.write(' \n')
            file.write("sp_geometry_xyz = wfn.molecule().to_string('xyz') \n")
            file.write(f"with open('{mol_name}_{job_name}.xyz', 'w') as xyz_file:\n")
            file.write("    xyz_file.write(sp_geometry_xyz) \n")
            file.write(' \n')
            file.write(f"psi4.core.set_active_molecule(wfn.molecule()) \n")
            file.write(' \n')
            file.write('orbital_energies = wfn.epsilon_a_subset("AO", "ALL").np \n')
            file.write(' \n')
            file.write('n_occ_alpha = wfn.nalpha() \n')
            file.write('homo_energy = orbital_energies[n_occ_alpha - 1] \n')
            file.write('lumo_energy = orbital_energies[n_occ_alpha] \n')
            file.write('energy_gap = lumo_energy - homo_energy \n')
            file.write(' \n')
            file.write('hartree_to_ev = 27.2114079527 \n')
            file.write('homo_energy_ev = homo_energy * hartree_to_ev \n')
            file.write('lumo_energy_ev = lumo_energy * hartree_to_ev \n')
            file.write('energy_gap_ev = lumo_energy_ev -  homo_energy_ev \n')
            file.write(' \n')
            file.write(f"with open('{mol_name}_{job_name}_energy_and_gap.txt', 'w') as file:\n")
            file.write('    file.write(f"Optimized energy: {energy:.6f} Hatree\\n") \n')
            file.write('    file.write(f"HOMO: {homo_energy_ev:.6f} eV\\n") \n')
            file.write('    file.write(f"LUMO: {lumo_energy_ev:.6f} eV\\n") \n')
            file.write('    file.write(f"Energy gap (HOMO-LUMO): {energy_gap_ev:.6f} eV\\n") \n')
        if (job_name=='opt'):
            file.write(f"energy, wfn = psi4.optimize('{functional}', return_wfn=True)  \n")
            file.write(' \n')
            file.write("optimized_geometry_xyz = wfn.molecule().to_string('xyz') \n")
            file.write(f"with open('{mol_name}_{job_name}.xyz', 'w') as xyz_file:\n")
            file.write("    xyz_file.write(optimized_geometry_xyz) \n")
            file.write(' \n')
            file.write(f"psi4.core.set_active_molecule(wfn.molecule()) \n")
            file.write(' \n')
            file.write('orbital_energies = wfn.epsilon_a_subset("AO", "ALL").np \n')
            file.write(' \n')
            file.write('n_occ_alpha = wfn.nalpha() \n')
            file.write('homo_energy = orbital_energies[n_occ_alpha - 1] \n')
            file.write('lumo_energy = orbital_energies[n_occ_alpha] \n')
            file.write('energy_gap = lumo_energy - homo_energy \n')
            file.write(' \n')
            file.write('hartree_to_ev = 27.2114079527 \n')
            file.write('homo_energy_ev = homo_energy * hartree_to_ev \n')
            file.write('lumo_energy_ev = lumo_energy * hartree_to_ev \n')
            file.write('energy_gap_ev = lumo_energy_ev -  homo_energy_ev \n')
            file.write(' \n')
            file.write(f"with open('{mol_name}_{job_name}_energy_and_gap.txt', 'w') as file:\n")
            file.write('    file.write(f"Optimized energy: {energy:.6f} Hatree\\n") \n')
            file.write('    file.write(f"HOMO: {homo_energy_ev:.6f} eV\\n") \n')
            file.write('    file.write(f"LUMO: {lumo_energy_ev:.6f} eV\\n") \n')
            file.write('    file.write(f"Energy gap (HOMO-LUMO): {energy_gap_ev:.6f} eV\\n") \n')
        file.write(' \n')


def write_psi4_reorg(job_name, mol_name, functional='b3lyp', basis_set='6-31g*'):
    """
    Generates a psi4 input file based on the provided parameters.
    
    Parameters
    ----------
    job_name (str): The type of job to run. Possible values:
        - 'cation': Geometry optimisation cation (opt_c) and single of neutral charge, cation geometry (n_c_geo)
        - 'anion': Geometry optimisation anion (opt_a) and single of neutral charge, anion geometry (n_a_geo)
        - 'sp_c': Single point calculation of neutral geometry at cation charge
        - 'sp_a': Single point calculation of neutral geometry at anion charge                                                                                                                                     
    mol_name (str): The name of the molecule.
    functional (str, optional): The functional to use. Default is 'b3lyp'.
    basis_set (str, optional): The basis set to use. Default is '6-31g'.
    
    Returns
    -------
         None: Writes the psi4 input file to disk.

    """
    if (job_name=='cation'):

        charge = 1
        mult = 2

        job_type = 'opt_c'

        geometry_file = f'{mol_name}_opt.xyz'
        geometry_sp_file = f'{mol_name}_opt_c.xyz'

        job_sp = 'n_c_geo'

    if (job_name=='anion'):

        charge = -1
        mult = 2

        job_type = 'opt_a'

        geometry_file = f'{mol_name}_opt.xyz'
        geometry_sp_file = f'{mol_name}_opt_a.xyz'

        job_sp = 'n_a_geo'

    if (job_name=='sp_c'):

        charge = 1
        mult = 2

        geometry_file = f'{mol_name}_opt.xyz'

    if (job_name=='sp_a'):

        charge = -1
        mult = 2

        geometry_file = f'{mol_name}_opt.xyz'

                                                                                          
    file_name = f'{mol_name}_{job_name}.py'
    #Includes information about basis set to allow for ramping
    memory_num = '4000 mb'

    with open(file_name, 'w') as file:
        file.write('import psi4\n')
        file.write(' \n')#
        file.write(f"psi4.set_memory('{memory_num}')\n")
        file.write('psi4.set_num_threads(4)\n')
        file.write(' \n')#
        file.write(f"mol = psi4.geometry(open('{geometry_file}').read())\n")#use geometry from opt
        file.write('psi4.core.set_active_molecule(mol)\n')#set as active molecule
        file.write(f"mol.set_molecular_charge({charge}) \n")# set charge
        file.write(f"mol.set_multiplicity({mult}) \n")# set multiplicity
        file.write(' \n')
        if (job_name=='cation') or (job_name=='anion'):
            file.write(f"psi4.set_options({{'basis': '{basis_set}', 'scf_type': 'uhf'}})\n")
            file.write(' \n')
            file.write(f"energy, wfn = psi4.optimize('{functional}', return_wfn=True)  \n")
            file.write(' \n')
            file.write("optimized_geometry_xyz = wfn.molecule().to_string('xyz') \n")
            file.write(f"with open('{mol_name}_{job_type}.xyz', 'w') as xyz_file:\n")
            file.write("    xyz_file.write(optimized_geometry_xyz) \n")
            file.write(' \n')
            file.write(f"psi4.core.set_active_molecule(wfn.molecule()) \n")
            file.write(' \n')
            file.write('orbital_energies = wfn.epsilon_a_subset("AO", "ALL").np \n')
            file.write(' \n')
            file.write('n_occ_alpha = wfn.nalpha() \n')
            file.write('homo_energy = orbital_energies[n_occ_alpha - 1] \n')
            file.write('lumo_energy = orbital_energies[n_occ_alpha] \n')
            file.write('energy_gap = lumo_energy - homo_energy \n')
            file.write(' \n')
            file.write('hartree_to_ev = 27.2114079527 \n')
            file.write('energy_ev =  energy * hartree_to_ev\n')
            file.write('homo_energy_ev = homo_energy * hartree_to_ev \n')
            file.write('lumo_energy_ev = lumo_energy * hartree_to_ev \n')
            file.write('energy_gap_ev = lumo_energy_ev -  homo_energy_ev \n')
            file.write(' \n')
            file.write(f"with open('{mol_name}_{job_type}_energy_and_gap.txt', 'w') as file:\n")
            file.write('    file.write(f"Optimized energy: {energy_ev:.6f} eV\\n") \n')
            file.write('    file.write(f"HOMO: {homo_energy_ev:.6f} eV\\n") \n')
            file.write('    file.write(f"LUMO: {lumo_energy_ev:.6f} eV\\n") \n')
            file.write('    file.write(f"Energy gap (HOMO-LUMO): {energy_gap_ev:.6f} eV\\n") \n')
            file.write(' \n')
            file.write(' \n') #n_c_geo or n_a_geo
            file.write(f"mol_sp = psi4.geometry(open('{geometry_sp_file}').read())\n")#use geometry from opt
            file.write('psi4.core.set_active_molecule(mol_sp)\n')#set as active molecule
            file.write(f"mol_sp.set_molecular_charge(0) \n")# set charge
            file.write(f"mol_sp.set_multiplicity(1) \n")# set multiplicity
            file.write(' \n')
            file.write(f"psi4.set_options({{'basis': '{basis_set}', 'scf_type': 'df'}})\n")
            file.write(' \n')
            file.write(f"energy, wfn = psi4.energy('{functional}', return_wfn=True)  \n")
            file.write(' \n')
            file.write(' \n')
            file.write(f"psi4.core.set_active_molecule(wfn.molecule()) \n")
            file.write(' \n')
            file.write('orbital_energies = wfn.epsilon_a_subset("AO", "ALL").np \n')
            file.write(' \n')
            file.write('n_occ_alpha = wfn.nalpha() \n')
            file.write('homo_energy = orbital_energies[n_occ_alpha - 1] \n')
            file.write('lumo_energy = orbital_energies[n_occ_alpha] \n')
            file.write('energy_gap = lumo_energy - homo_energy \n')
            file.write(' \n')
            file.write('energy_ev =  energy * hartree_to_ev\n')
            file.write('homo_energy_ev = homo_energy * hartree_to_ev \n')
            file.write('lumo_energy_ev = lumo_energy * hartree_to_ev \n')
            file.write('energy_gap_ev = lumo_energy_ev -  homo_energy_ev \n')
            file.write(' \n')
            file.write(f"with open('{mol_name}_{job_sp}_energy_and_gap.txt', 'w') as file:\n")
            file.write('    file.write(f"Single Point energy: {energy_ev:.6f} eV\\n") \n')
            file.write('    file.write(f"HOMO: {homo_energy_ev:.6f} eV\\n") \n')
            file.write('    file.write(f"LUMO: {lumo_energy_ev:.6f} eV\\n") \n')
            file.write('    file.write(f"Energy gap (HOMO-LUMO): {energy_gap_ev:.6f} eV\\n") \n')
            file.write(' \n')
        if (job_name=='sp_c') or (job_name=='sp_a'):
            file.write(f"psi4.set_options({{'basis': '{basis_set}', 'scf_type': 'uhf'}})\n")
            file.write(' \n')
            file.write(f"energy, wfn = psi4.energy('{functional}', return_wfn=True)  \n")
            file.write(' \n')
            file.write(' \n')
            file.write(f"psi4.core.set_active_molecule(wfn.molecule()) \n")
            file.write(' \n')
            file.write('orbital_energies = wfn.epsilon_a_subset("AO", "ALL").np \n')
            file.write(' \n')
            file.write('n_occ_alpha = wfn.nalpha() \n')
            file.write('homo_energy = orbital_energies[n_occ_alpha - 1] \n')
            file.write('lumo_energy = orbital_energies[n_occ_alpha] \n')
            file.write('energy_gap = lumo_energy - homo_energy \n')
            file.write(' \n')
            file.write('hartree_to_ev = 27.2114079527 \n')
            file.write('energy_ev =  energy * hartree_to_ev\n')
            file.write('homo_energy_ev = homo_energy * hartree_to_ev \n')
            file.write('lumo_energy_ev = lumo_energy * hartree_to_ev \n')
            file.write('energy_gap_ev = lumo_energy_ev -  homo_energy_ev \n')
            file.write(' \n')
            file.write(f"with open('{mol_name}_{job_name}_energy_and_gap.txt', 'w') as file:\n")
            file.write('    file.write(f"Single Point energy: {energy_ev:.6f} eV\\n") \n')
            file.write('    file.write(f"HOMO: {homo_energy_ev:.6f} eV\\n") \n')
            file.write('    file.write(f"LUMO: {lumo_energy_ev:.6f} eV\\n") \n')
            file.write('    file.write(f"Energy gap (HOMO-LUMO): {energy_gap_ev:.6f} eV\\n") \n')
            file.write(' \n')

            