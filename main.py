from fragments_io import *
from torsion_io import *
from write_input_io import *

def submit_slurm_job(dimer_name):
    """
    Submit a slurm job

    """

    string = f'sbatch {dimer_name}.sh'

    process = subprocess.run(string,
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE, shell=True,check=True)
    return process.stdout

def torsional_scan_run(example_mols_pres):

    #Loads the smi file and creates a dictionary of fragments
    farg_dic = fragment_dictionary('example_mols_pres')
    #Creates a dictionary of all the dimers from the fragment dictionary
    dimers_dic = combined_fragments_dic(farg_dic)

    #creates a directory named tortional_scan
    os.mkdir('tortional_inputs')

    #goes into that directory
    os.chdir('tortional_inputs')

    #torsion_dic, bond_dic, mol_dic = {}
    for dimer_name, dimer_smiles in dimers_dic.items():
    #makes directory called after the job name
        os.mkdir(f'{dimer_name}-job')
        #goes into that directory
        os.chdir(f'{dimer_name}-job')

        #turns smiles string into rdkit object
        mol = Chem.MolFromSmiles(dimer_smiles)
        #finds the bond between the fragment
        bond = getBond(mol)
        #torsion of the bond between the fragment
        torsion = getTorsion(mol, bond[0])
        #gets rdkit estimated coordinates of dimer
        mold3d = embed_molecule(mol)

        #writes a guassian input file
        write_gaussian_scan_input(mold3d, dimer_name, torsion)
        #writes the slurm file
        write_slurm_input(dimer_name)
        #submits the slurm jon
        submit_slurm_job(dimer_name)
        #goes back to previous directory
        os.chdir(os.path.dirname(os.getcwd()))
