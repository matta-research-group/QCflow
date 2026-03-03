import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import numpy as np
import pandas as pd
from QCflow.load_gaussian import *
from QCflow.torsion_parser import *
from QCflow.find_torsion import *
from QCflow.write_psi4 import *
from QCflow.run_psi4 import *
from QCflow.energy_calculations import *
import os

# This is an example for one molecule
# This can easily be converted into a loop for multiple molecules
# The getTorsion assumes single bonded but can be altered if required
# If you have more than 1 torsion or tiple, double etc bond link
# Refer to this --> https://github.com/matta-research-group/QCflow/blob/qcflow-0.5/QCflow/torsion_parser.py
# torsion_parser.py also has the scripts to extract the data you may want
# The docs  --> https://matta-research-group.github.io/QCflow/

smiles_example = 'CC(=O)Oc1c(OC(C)=O)c(-c2c(OC(C)=O)c(OC(C)=O)cc3cc[nH]c23)c2cc[nH]c2c1'

mol = Chem.MolFromSmiles(smiles_example)

#adds H's and 3D coordinates
mol = embed_molecule(mol)

#finds the rotatable bond in the molecule
bond = getBond(mol)

#Gets the whole torsion
mol_torsions = getTorsion(mol, bond[0])

#writes the Gaussian input file for the torsion scan
write_gaussian('tor', 'example_mol', smiles_example, functional='B3LYP', basis_set='6-31G*', mol=mol, 
                                                                                       torsion=mol_torsions)
#writes the SLURM script to run the Gaussian job
write_slurm('tor', 'example_mol', cpus=10)

#submits the SLURM job, with retries in case of failure
submit_slurm_job('tor', 'example_mol', max_retries=5, wait_seconds=30)