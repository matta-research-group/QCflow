import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import numpy as np
import pandas as pd
import QCflow
from QCflow.load_gaussian import *
from QCflow.torsion_parser import *
from QCflow.find_torsion import *
from QCflow.write_psi4 import *
from QCflow.run_psi4 import *
from QCflow.energy_calculations import *
import os

df = pd.read_csv('FilteredSmilesFragments_props.csv')

smiles_dict = {str(i): row['smiles'] for i, row in df.iterrows()}

for k, v in smiles_dict.items():
    run_psi4('opt', k, v, 12, 10, 'b3lyp', '3-21G')
