from fragments import *
from torsion_run import *
from write_gaussian import *
from load_gaussian import *
from slurm import *
import json
import csv
import matplotlib.pyplot as plt
import cclib

def save_torsion(data, mol_name):
    """
    Saves the torsional profile of a dimer

    data : Loaded data from the .log file containing the torsional information

    mol_name : number name of dimer or trimer
    """
    data.scanenergies = data.scfenergies[data.optstatus == 4]

    x = data.scanparm

    y = data.scanenergies-np.min(data.scanenergies)

    np.savetxt(f'{mol_name}_tor_prof.csv', np.vstack((x,y)).T, delimiter=', ')

def find_min_energy(data):
    """
    Finds the miniumum energy of the torsion

    data : Loaded data from the .log file containing the torsional information
    """
    #finds opt energy of each 10 deg scan
    data.scanenergies = data.scfenergies[data.optstatus == 4]
    #finds min energy of all opt 36 scans
    mini = np.where(data.scanenergies==np.min(data.scanenergies))
    #turns location of the miniumum energy into a number
    mini_value = int(mini[-1])

    return mini_value

def torsional_parser(mol_name, mol_dic):
    """
    When provided with the name of the molecule and the molecule dictionary
    that it came from this function saves the torsional profile as a .csv file.
    In addition it also returns the updated geometry and the dihedral of the
    miniumum

    mol_name : number name of dimer or trimer

    mol_dic : dictionary in which the dimer or trimer is in
    """
    mol = Chem.MolFromSmiles(mol_dic[mol_name])
    #Converts the SMILE into rdkit readable string
    mol_h = AllChem.AddHs(mol)
    #Adds hydrogen atoms
    mol_3d = embed_molecule(mol_h)
    #Embeds the molecule
    conf = mol_3d.GetConformer()
    #Turns into a conformer
    data = load_torsional_data(mol_name)
    #Loads the .log file containing the torsional data
    save_torsion(data, mol_name)
    #Saves the torsional scan as .csv file with the name of the mol
    min_energy = find_min_energy(data)
    #Finds the minimum energy torsion
    min_ang_loc = np.where(data.scanenergies==np.min(data.scanenergies))
    #Location of the min angle
    min_angle = data.scanparm[0][min_ang_loc[0][0]]
    #Dihedral angle of the min energy

    for i in range(conf.GetNumAtoms()):
        correct_pos = conf.SetAtomPosition(i, data.converged_geometries[min_energy][i])
    #Uses this minium energy torsion to give the correct geometry
    return conf, min_angle
