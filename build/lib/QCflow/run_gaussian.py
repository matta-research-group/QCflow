from QCflow.fragments import *
from QCflow.find_torsion import *
from QCflow.write_gaussian import *
from QCflow.slurm import *
from QCflow.torsion_parser import *
import os
import json

def rdkit_predict_conf(mol_smiles, num_of_conformer=100, max_iter=500, min_energy_MMFF=10000, min_energy_index_MMFF=0):
    """
    Generates conformers for a given molecule using RDKit and returns the lowest energy conformer.

    Parameters:
    mol_smiles (str): The SMILES representation of the molecule.
    num_of_conformer (int): The number of conformers to be generated (default: 100).
    max_iter (int): The maximum number of iterations for conformer optimization (default: 500).
    min_energy_MMFF (float): The minimum energy threshold for selecting the lowest energy conformer (default: 10000).
    min_energy_index_MMFF (int): The index of the lowest energy conformer (default: 0).

    Returns:
    conf (Chem.Conformer): The lowest energy conformer of the molecule.

    """
    # Number of conformers to be generated
    #num_of_conformer=100
    #max_iter=500
    # Default values for min energy conformer
    #min_energy_MMFF=10000
    #min_energy_index_MMFF=0
    
    mol = Chem.MolFromSmiles(mol_smiles)
    mol_h_MMFF = Chem.AddHs(mol)
    cids = AllChem.EmbedMultipleConfs(mol_h_MMFF, numConfs=num_of_conformer,params=AllChem.ETKDG()) # Generate conformers (stored in side the mol object)
    results_MMFF = AllChem.MMFFOptimizeMoleculeConfs(mol_h_MMFF,maxIters=max_iter)

    for index, result in enumerate(results_MMFF):
        if(min_energy_MMFF>result[1]):       
            min_energy_MMFF=result[1]
            min_energy_index_MMFF=index

    lowest_mol = Chem.Mol(mol_h_MMFF,False,min_energy_index_MMFF)

    conf = lowest_mol.GetConformer()

    return conf


def run_calc(job_name, mol_name, mol_smile, functional='B3LYP', basis_set='6-31G*'):
    """
    Submits a gaussian calculation to CREATE HPC.

    Parameters:
    job_name (str): Name of the job to be submitted.
    mol_name (str): Name of the oligomer.
    mol_smile (str): SMILE string of the oligomer.
    functional (str, optional): Quantum chemistry functional to be used (default is 'B3LYP').
    basis_set (str, optional): Basis set to be used (default is '6-31G*').

    This function performs the following steps:
    1. Creates a directory for the molecule.
    2. Changes the current working directory to the newly created directory.
    3. Converts the SMILES string to an RDKit molecule object.
    4. Generates 3D coordinates for the molecule.
    5. Writes a Gaussian input file for the molecule.
    6. Writes a SLURM script for job submission.
    7. Submits the SLURM job.
    8. Returns to the previous directory.
    """

    #makes a directory for the molecule
    os.mkdir(f'{mol_name}')
    #goes into directory
    os.chdir(f'{mol_name}')
    #turns smiles string into rdkit object
    mol = Chem.MolFromSmiles(mol_smile)
    #gets rdkit estimated coordinates of dimer
    mol3d = embed_molecule(mol)

    conf_geo = rdkit_predict_conf(mol_smile)
    #writes a guassian input file
    write_gaussian(job_name, mol_name, mol_smile, functional, basis_set, mol=mol3d, torsion=0, conformer=conf_geo)
    #writes the slurm file
    write_slurm(job_name, mol_name)
    #submits the slurm jon
    submit_slurm_job(job_name, mol_name)
    #goes back to previous dirctory
    os.chdir(os.path.dirname(os.getcwd()))

def staging_opt(job_name, mol_name, mol_smile, mol_dic, functional, basis_set):
    """
    Checks if the calcultion has been completed at has been a success.
    If it has, then appends a dictionary showing this. If the calculations haven't been run
    all the way, then runs them. If the calculation has failed, then appends a dictionary
    
    Parameters:
    job_name (str): Type of job run e.g. pop_opt_n
    mol_name (str): The name of the oligomer as seen in the dictionary i.e. if melanin fragment (b) is combined
    mol_smile (str): SMILE string of oligomer
    mol_dic (dict): Dictionary of oligomers where key is the name of the oligomer and value is the SMILES string
    functional (str): Functional used in calculations (e.g. B3LYP)
    basis_set (str): Basis set used (e.g. 6-31G*)
    
    Returns:
    tuple: A tuple containing three dictionaries:
        fully_complete (dict): Oligomers that have been calculated at the highest basis set
        not_complete (dict): Oligomers that failed and need manual assessment (shows basis set they failed at)
        in_progress (dict): Oligomers that are still in progress (shows basis set they are currently being run at)
    """

    fully_complete = {}
    not_complete = {}
    in_progress = {}
    k = mol_name

    if os.path.isfile(f'{k}/{k}_{job_name}.log') == False:
        #Then run it
        run_calc(job_name, k, mol_smile, mol_dic, functional, basis_set)
        in_progress[f'{k}'] = basis_set

    if os.path.isfile(f'{k}/{k}_{job_name}.log') == True:
    #If the file exists
        data = cclib.ccopen(f'{k}/{k}_{job_name}.log').parse()
        #parse that file
        if data.metadata['success'] == False:
            not_complete[f'{k}'] = basis_set
            #add to fail dictionary
        if data.metadata['success'] == True:
            #If the calculations worked, was it done at correct basis set
            if data.metadata['basis_set'] == basis_set:
        
                fully_complete[f'{k}'] = basis_set

            if data.metadata['basis_set'] != basis_set:
                #Run the calc at correct basis set
                run_calc(job_name, k, mol_smile, mol_dic, functional, basis_set)
                in_progress[f'{k}'] = basis_set


    return fully_complete, not_complete, in_progress


def run_torsion(mol_name, mol_smiles, functional='B3LYP', basis_set='6-31G*'):
    """
    When provided with a dictionary of molecules, this function creates the SLURM 
    and Gaussian input files and then submits the torsional scan.

    Parameters:
    mol_name (str): Name of the oligomer.
    mol_smile (str): SMILE string of the oligomer.
    functional (str): The functional to be used in Gaussian calculations. Default is 'B3LYP'.
    basis_set (str): The basis set to be used in Gaussian calculations. Default is '6-31G*'.

    The function performs the following steps for each molecule:
    1. Creates a directory named after the molecule.
    2. Converts the SMILES string to an RDKit molecule object.
    3. Identifies the bond and torsion of the molecule.
    4. Generates 3D coordinates for the molecule.
    5. Writes the Gaussian input file.
    6. Writes the SLURM job file.
    7. Submits the SLURM job.
    8. Returns to the previous directory.
    """

    
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