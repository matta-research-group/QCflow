from QCflow.fragments import *
from QCflow.find_torsion import *
from QCflow.write_gaussian import *
from QCflow.slurm import *
from QCflow.torsion_parser import *
import os
import json


def run_opt_neutral(mol_name, mol_smile, mol_dic, functional='B3LYP', basis_set='6-31G*'):
    """
    Submits a neutral optimization with population analysis when provided with the name of the mol

    mol_name : name of oligomer

    mol_smile : SMILE string of oligomer

    mol_dic : Dictionary the oligomer came from

    functional : preset is B3LYP

    basis_set : preset is 6-31G*
    """
    #goes into directory
    os.chdir(f'{mol_name}')
    #turns smiles string into rdkit object
    mol = Chem.MolFromSmiles(mol_smile)
    #gets rdkit estimated coordinates of dimer
    mol3d = embed_molecule(mol)
    #Saves the torsional scan as .csv and finds the lowest energy geometry
    conf_geo = torsion_parser(mol_name, mol_dic)
    #writes a guassian input file
    write_gaussian('pop_opt_n', mol_name, mol_smile, functional, basis_set, mol3d, 0, conf_geo)
    #writes the slurm file
    write_slurm('pop_opt_n', mol_name)
    #submits the slurm jon
    submit_slurm_job('pop_opt_n', mol_name)
    #goes back to previous directory
    os.chdir(os.path.dirname(os.getcwd()))


def staging_opt(mol_name, mol_smile, mol_dic, job_type, functional, basis_set1, basis_set2):
    """
    Checks if the geometry optimization as been completed at the highest desried basis set.
    If it has then appends a dictionary showing this. If the calculations haven't been ran
    all the way then runs them. If the calculation has failed then appends a dictionary
    with the name of the oligomer and the basis set it failed at.

    mol_name : the name of the oligomer as seen in the dictionary i.e. if melanin fragment (b) is combined
                with organic electronic fragment (1) in the v1 position will be b_1_v1

    mol_smile : SMILE string of oligomer

    mol_dic : dictionary of oligomers where key is the name of the oligomer and value is the SMILES string

    job_type : Type of job ran e.g. pop_opt_n

    functional : Functional used in calculations (e.g. B3LYP)

    basis_set1 : Lowest basis set used (e.g.  3-21G*)

    basis_set2 : Highest basis set used Preset (e.g. 6-31G*)

    Returns three dictionaries:

        fully_complete : Oligomers that have been calculated at the highest basis set

        not_complete : Oligomers that failed and need manual assessment (shows basis set they failed at)

        in_progress : Oligomers that are still in progress (shows basis set they are currently being ran at)


    """
    fully_complete = {}
    not_complete = {}
    in_progress = {}
    k = mol_name

    if os.path.isfile(f'{k}/{k}_{job_type}.log') == False:
        #Then run it
        run_opt_neutral(k, mol_smile, mol_dic, functional, basis_set1)
        in_progress[f'{k}'] = basis_set1

    if os.path.isfile(f'{k}/{k}_{job_type}.log') == True:
    #If the file exists
        data = cclib.ccopen(f'{k}/{k}_{job_type}.log').parse()
        #parse that file
        if data.metadata['success'] == False:
            not_complete[f'{k}'] = basis_set1
            #add to fail fictionary
        if data.metadata['success'] == True:
            #If the calculations worked, was it done at highest basis set
            if data.metadata['basis_set'] == basis_set1:
                #Run the calc at higher basis set if no
                run_opt_neutral(k, mol_smile, mol_dic, functional, basis_set2)
                in_progress[f'{k}'] = basis_set2

            if data.metadata['basis_set'] == basis_set2:
                #If ran at highest basis set then calculation complete
                fully_complete[f'{k}'] = basis_set2


    return fully_complete, not_complete, in_progress
