from fragments_io import *
from torsion_io import *
from write_input_io import *
from file_io import *


def run_torsion(fragment_smi_file):

    '''
    When provided with a .smi file with all the fragments, this function
    combines all the fragemnts into dimers and then creates the SLURM and
    guassian input file and then submits the jobs

    fragment_smi_file : name of the .smi file containing the fragments
    '''

    #Loads the smi file and creates a dictionary of fragments
    frag_dic = read_fragments('fragment_smi_file')
    #Creates a dictionary of all the dimers from the fragment dictionary
    mol_dic = make_dimers_dic(frag_dic)

    #creates a directory named tortional_scan
    os.mkdir('torsional_inputs')

    #goes into that directory
    os.chdir('torsional_inputs')

    for mol_name, mol_smiles in mol_dic.items():
    #makes directory called after the job name
        os.mkdir(f'{mol_name}-job')
        #goes into that directory
        os.chdir(f'{mol_name}-job')

        #turns smiles string into rdkit object
        mol = Chem.MolFromSmiles(mol_smiles)
        #finds the bond between the fragment
        bond = getBond(mol)
        #torsion of the bond between the fragment
        torsion = getTorsion(mol, bond[0])
        #gets rdkit estimated coordinates of dimer
        mol3d = embed_molecule(mol)

        #writes a guassian input file
        write_gaussian('torsional scan', mol_name, mol3d, torsion)
        #writes the slurm file
        write_slurm('torsional scan', mol_name)
        #submits the slurm jon
        submit_slurm_job('torsional scan', mol_name)
        #goes back to previous directory
        os.chdir(os.path.dirname(os.getcwd()))
