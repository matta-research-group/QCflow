from fragments import *
from find_torsion import *
from write_gaussian import *
from slurm import *
from torsion_parser import *


def run_torsion(fragment_smi_file, functional='B3LYP', basis_set='6-31G*'):

    '''
    When provided with a .smi file with all the fragments, this function
    combines all the fragemnts into dimers and then creates the SLURM and
    guassian input file and then submits the jobs

    fragment_smi_file : name of the .smi file containing the fragments

    functional : Preset is B3LYP

    basis_set : Preset is 6-31G*
    '''

    #Loads the smi file and creates a dictionary of fragments
    frag_dic = read_fragments(fragment_smi_file)
    #Creates a dictionary of all the dimers from the fragment dictionary
    mol_dic = make_dimer_dic(frag_dic)
    #saves the dictionary as a json file
    save_dictionary(mol_dic, 'mol_dic.json')

    for mol_name, mol_smiles in mol_dic.items():
    #makes directory called after the job name
        os.mkdir(f'{mol_name}')
        #goes into that directory
        os.chdir(f'{mol_name}')

        #turns smiles string into rdkit object
        mol = Chem.MolFromSmiles(mol_smiles)
        #finds the bond between the fragment
        bond = getBond(mol)
        #torsion of the bond between the fragment
        torsion = getTorsion(mol, bond[0])
        #gets rdkit estimated coordinates of dimer
        mol3d = embed_molecule(mol)

        #writes a guassian input file
        write_gaussian('tor', mol_name, mol_smiles, functional, basis_set, mol3d, torsion)
        #writes the slurm file
        write_slurm('tor', mol_name)
        #submits the slurm jon
        submit_slurm_job('tor', mol_name)
        #goes back to previous directory
        os.chdir(os.path.dirname(os.getcwd()))
