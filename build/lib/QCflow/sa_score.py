import sys
import os
import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(True)
sys.path.append(os.path.join(os.environ['CONDA_PREFIX'],'share','RDKit','Contrib'))
from SA_Score import sascorer

def sa_scorer(smile):
    """
    Carries out a synthetic accessibility score on a molecule (SA score)
    The literature for this paper can be found here: https://doi.org/10.1186/1758-2946-1-8

    Parameters
    ----------
    smile (str): The SMILES string of the molecule to be scored

    Returns
    -------
    sa_score_val (float): The SA score of the molecule
    """
    m = Chem.MolFromSmiles(smile)

    #Run synethic accessibility score
    sa_score_val = sascorer.calculateScore(m)
    return sa_score_val