import rdkit
import itertools
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import itertools
from itertools import combinations
from typing import List
from CombineMols.CombineMols import CombineMols

def adding_attach(smi, find='[cH;x2]', get_rid='C([*])'):
    """
    Adds an attachment point to a given fragment.

    Parameters:
    smi (str): The SMILES string of the fragment.
    find (str): The SMARTS pattern of the possible attachment point. Default is '[cH;x2]'.
    get_rid (str): The replacement pattern for the attachment point. Default is 'C([*])'.

    Returns:
    list: A list of unique SMILES strings representing the fragment with the added attachment point.
    """
    finder = Chem.MolFromSmarts(find)
    replace = Chem.MolFromSmiles(get_rid)
    mol = Chem.MolFromSmiles(smi)

    one_attach = rdkit.Chem.rdmolops.ReplaceSubstructs(mol, finder, replace)

    one_attach_smi = []
    for x in one_attach:
        smi = Chem.MolToSmiles(x)
        can_smi = Chem.CanonSmiles(smi)
        one_attach_smi.append(can_smi)

    cleaned = list(set(one_attach_smi))

    return cleaned

def generate_attachment_points(frag_dic):
    """
    This function takes a dictionary of fragments and generates attachment points for each fragment.
    The attachment points are determined by using the 'adding_attach' function.
    The attachment points are represented as alphabetical characters ('A', 'B', ...) and are appended to the identifiers of the fragments.
    The resulting attachment points dictionary is returned.

    Parameters:
    frag_dic (dict): A dictionary containing fragment information with keys as identifiers and values as SMILE strings.

    Returns:
    dict: A dictionary containing attachment points. Keys are composed of identifiers followed by
          alphabetical characters ('A', 'B', ...) representing different attachment points for each fragment.
          Values are the corresponding attachment points obtained from the 'adding_attach' function. They are
          a rdkit mol object.

    """

    attach_dic = {}

    p = Chem.MolFromSmiles('I')

    for k, v in frag_dic.items():
        all_attach = adding_attach(v, find='[cH;^2]', get_rid='C([I])') #using 'I' because a combine function uses it later
        #finding sp2 carbons
        mol_test = Chem.MolFromSmiles(all_attach[0])

        if mol_test.HasSubstructMatch(p) == True:

            for i, attach_point in enumerate(all_attach):
                attach_dic[f'{k}_{chr(ord("A") + i)}'] = attach_point #adds letter as some have more than one attachment

        if mol_test.HasSubstructMatch(p) == False: #if there is no iodine attachment, add hydrogens and try aliphatic carbons
            add_h = Chem.AddHs(mol_test)
            smi_h = Chem.MolToSmiles(add_h)
            all_attach2 = adding_attach(smi_h, find='[CH;^2]', get_rid='C([I])') #attaching at aliphatic carbons

            for i, attach_point in enumerate(all_attach2):
                attach_dic[f'{k}_{chr(ord("A") + i)}'] = attach_point #adds letter as some have more than one attachment

    return attach_dic

def combine_structure(molecule_a, molecule_b):
    """
    Combines two molecules together at a specified attachment point.

    Parameters:
    molecule_a (object): The first molecule to be combined.
    molecule_b (object): The second molecule to be combined.

    Returns:
    object: The combined molecule.

    """
    new_mol = CombineMols(molecule_a, molecule_b, "I")
    # Combining two dimers together at attachment point denoted atom 'I'

    new_mol_final = new_mol[0]

    new_mol_smi = Chem.MolToSmiles(new_mol_final)

    return new_mol_smi

def make_molecule_dic_from_2_dic(fragment_dic_1, fragment_dic_2):
    """
    Generates a dictionary of molecule by combining two different fragment dictionaries.

    Parameters:
        fragment_dic_1 (dict): The first fragment dictionary.
        fragment_dic_2 (dict): The second fragment dictionary.

    Returns:
        dict: A dictionary of dimers, where the key is the name of the molecule and the value is the SMILES string.

    """

    mol_smiles = [combine_structure(v1, v2)
                    for v1, v2 in itertools.product(list(fragment_dic_1.values()),list(fragment_dic_2.values()))]

    mol_names = [(f'{k1}_{k2}')
                    for k1, k2 in itertools.product(fragment_dic_1.keys(),fragment_dic_2.keys())]

    mol_dic = { k : v for k, v in zip(mol_names, mol_smiles) }

    return mol_dic