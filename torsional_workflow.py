import numpy as np
import matplotlib.pyplot as plt
import rdkit
import cclib
import itertools
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import requests
from itertools import combinations
from typing import List
import os
import os.path
import shutil
import subprocess

def combine_fragments(fragments: List[str]) -> List[Chem.rdchem.Mol]:
    """
    Gets all SMILES for all possible combinations (of 2 fragments) of the fragment list
    Creates a list of all the possible combinations as RDKit info
    """

    all_smiles = (f'{smiles1}.{smiles2}'.replace('(*)', '9')
                    for smiles1, smiles2 in combinations(fragments, r=2))
    # Convert all SMILES to mol objects
    mols = [Chem.MolFromSmiles(smi) for smi in all_smiles]
    return mols

def dimer_product(list1, list2):
    """ combines any fragments from 2 lists into dimers
    return the corresponding rdkit molecule
    """
    dimer_smiles = (f'{smiles1}.{smiles2}'.replace('(*)', '9')
                  for smiles1, smiles2 in itertools.product(list1,list2))

    dimers = [Chem.MolFromSmiles(smi) for smi in dimer_smiles]
    return dimers

def getBond(mol):
    """
    return rotatable bond using SMARTS pattern
    todo: fix for exceptions and reimplement
    """
    pattern = Chem.MolFromSmarts('[R!$(*#*)&!D1]-!@[R!$(*#*)&!D1]')
    bonds = mol.GetSubstructMatches(pattern)
    return bonds

#3 14 21 22 S 35 10.0
def getTorsion(mol,bond):
    """
    given a bond, returns a dihedral to scan
    """
    # get neighbors of first atom in bond
    for atom in mol.GetAtomWithIdx(bond[0]).GetNeighbors():
        idx = atom.GetIdx()
        mass = atom.GetMass()
        if idx!=bond[1]: #excludes the other atom in the bond
            if mass > 13: # N, S, O get priority
                first=idx
                break
            if mass > 12: # otherwise C
                first=idx
    # get neighbors of second atom in bond
    for atom in mol.GetAtomWithIdx(bond[1]).GetNeighbors():
        idx=atom.GetIdx()
        mass = atom.GetMass()
        if idx!=bond[0]: #excludes the other atom in the bond
            if mass > 13: # N, S, O get priority
                last=idx
                break
            if mass > 12: # otherwise C
                last=idx

    return (first, bond[0], bond[1], last)

def embed_molecule(mol):
    """
    takes a 2D rdkit molecule
    returns a 3D embedded rdkit molecule
    """
    addhs = Chem.AddHs(mol)
    AllChem.EmbedMolecule(addhs) # the embedded molecule
    return addhs

def write_gaussian_scan_input(dimer, dimer_name, dihedral_atoms, npoints=35, step=10):
    """
    writes gaussian input file for torsional scan
    default method is B3LYP/6-31G*
    npoints = number of optimisation points along the PES
    step = degrees between each scan point
    """
    # get atom names
    symbols = [a.GetSymbol() for a in dimer.GetAtoms()]
    # get x y z coords
    geometry=dimer.GetConformers()[0]
    # unpack torsion terms and add 1 for fortran 1-based numbering
    a1, a2, a3, a4 = [t+1 for t in dihedral_atoms]
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
        file.write(f'{a1} {a2} {a3} {a4} S {npoints} {nstep}.0\n')
        file.write(' \n')

def write_slurm_input(dimer_name):
    """
    creates slurm submission script named dimer_name.sh
    """

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
        file.write(f'module load test_switch_kcl/1.0.0-gcc-9.4.08. \n')
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

def submit_slurm_job(dimer_name):
    """
    Submit a slurm job
    """

    string = f'sbatch {dimer_name}.sh'

    process = subprocess.run(string,
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE, shell=True,check=True)
    return process.stdout

def delete_sh(dimer_name):
    return os.remove(f'{dimer_name}.sh')

def delete_com(dimer_name):
    return os.remove(f'{dimer_name}.com')

def TortionalScanAuto(smiles_list_with_attach):

    dimers_list = combine_fragments(smiles_list_with_attach) #all the possible combinations of fragments

    bonds_list = []
    for dimer in dimers_list:
        bonds_list.append(getBond(dimer)) #finds the bonds between the fragments

    torsions_list=[]
    for bond, dimer in zip(bonds_list, dimers_list):
        torsions_list.append(getTorsion(dimer,list(bond[0]))) #obtains the torsion of the bond between fragments

    dimers_3dlist=[]
    for dimer in dimers_list:
        dimer = embed_molecule(dimer)
        dimers_3dlist.append(dimer) #get coordinates of the dimer

    dimer_dic = { i : d for i, d in enumerate(dimers_3dlist)} #creates a list and numbers all the dimers

    os.mkdir('tortional_inputs') #creates a directory named tortional_scan

    os.chdir('tortional_inputs') #goes into that directory

    for (k, v), t in zip(dimer_dic.items(), torsions_list):
        os.mkdir(f'{k}-Job')
        os.chdir(f'{k}-Job')
        write_gaussian_scan_input(v, k, t) #writes the input file
        write_slurm_input(k) #writes the SLURM File
        submit_slurm_job(k) #submits SLURM jobs
        os.chdir(os.path.dirname(os.getcwd())) #goes back to previous directory
        #delete_sh(k) #deletes all the SLURM files
        #delete_com(k) #deletes all the inputs
