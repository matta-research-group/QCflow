from QCflow.fragments import *
from QCflow.torsion_run import *
from QCflow.write_gaussian import *
from QCflow.load_gaussian import *
from QCflow.slurm import *
import json
import csv
import matplotlib.pyplot as plt
import cclib

def save_torsion(data, mol_name):
    """
    Saves the torsional profile of a dimer

    data : Loaded cclib data from the .log file containing the torsional information

    mol_name : number name of dimer or trimer
    """
    data.scanenergies = data.scfenergies[data.optstatus == 4]

    x = data.scanparm

    y = data.scanenergies-np.min(data.scanenergies)

    np.savetxt(f'{mol_name}_tor_prof.csv', np.vstack((x,y)).T, delimiter=', ')

def find_min_energy(data):
    """
    Finds the minimum energy of the torsion

    data : Loaded cclib data from the .log file containing the torsional information
    """
    #finds opt energy of each 10 deg scan
    data.scanenergies = data.scfenergies[data.optstatus == 4]
    #finds min energy of all opt 36 scans
    mini = np.where(data.scanenergies==np.min(data.scanenergies))
    #gets index corresponding to minimum energy
    mini_value = int(mini[-1])

    return mini_value

def min_angle(data):
    """
    Finds the angle of the minimum energy torsion

    data : Loaded cclib data from the .log file containing the torsional information
    """
    data.scanenergies = data.scfenergies[data.optstatus == 4]
    min_ang_loc = np.where(data.scanenergies==np.min(data.scanenergies))
    #Location of the min angle
    min_ang = data.scanparm[0][min_ang_loc[0][0]]
    #Dihedral angle of the min energy
    return min_ang

def min_angle_dic(mol_name, min_ang):
    """
    Creates a dictionary of molecules with the angle of the minimum energy torsion

    mol_name : number name of dimer or trimer

    min_ang : angle of the minimum energy torsion in degrees
    """
    angle_dic = { k : v for k, v in zip(mol_name, min_ang) }

    return angle_dic

def save_angle(angle_dic, name_of_dic):
    """
    Saves the dictionary of all the minimum energy torsion

    angle_dic : dictionary of molecules with the angle of the minimum energy torsion

    name_of_dic : desried name of the json file i.e. name_of_dic.json
    """

    save_dictionary(angle_dic, name_of_dic)

def torsion_parser(mol_name, mol_dic):
    """
    When provided with the name of the molecule and the molecule dictionary
    that it came from this function saves the torsional profile as a .csv file.

    mol_name : number name of dimer or trimer
    mol_dic : dictionary in which the dimer or trimer is in

    Returns:
    optimised geometry at lowest energy minimum

    """
    mol = Chem.MolFromSmiles(mol_dic[mol_name])
    #Converts the SMILE into rdkit readable string
    mol_h = AllChem.AddHs(mol)
    #Adds hydrogen atoms
    mol_3d = embed_molecule(mol_h)
    #Embeds the molecule
    conf = mol_3d.GetConformer()
    #Turns into a conformer
    data = load_data(mol_name, 'tor')
    #Loads the .log file containing the torsional data
    save_torsion(data, mol_name)
    #Saves the torsional scan as .csv file with the name of the mol
    min_energy = find_min_energy(data)
    #Finds the minimum energy torsion

    for i in range(conf.GetNumAtoms()):
        correct_pos = conf.SetAtomPosition(i, data.converged_geometries[min_energy][i])
    #Uses this min energy torsion to give the correct geometry
    return conf

def setting_dihedral(mol_smile, deg):
    """
    Settings the dihedral angle of the torsion for individual scans

    mol_smiles : smi string of dimer

    deg : desired dihedral angle (0 or 180 for planar)

    returns a dimer with a dihedral angle that is inputted by user
    """
    mol = Chem.MolFromSmiles(mol_smile)
    #turns smiles string into rdkit object
    mol3d = embed_molecule(mol)
    #gets rdkit estimated coordinates of dimer
    conf = mol3d.GetConformer()
    #getting the conformer

    bond = getBond(mol)
    #finds the bond between the fragment
    torsion = getTorsion(mol, bond[0])
    #torsion of the bond between the fragment
    t1, t2, t3, t4 = [t for t in torsion]
    #getting atoms involved in torsion
    rdkit.Chem.rdMolTransforms.SetDihedralDeg(conf, t1, t2, t3, t4, deg)
    #setting the dihedral angle
    mp = AllChem.MMFFGetMoleculeProperties(mol3d)
    #mol properties
    ff = AllChem.MMFFGetMoleculeForceField(mol3d, mp)
    #force field
    for i in torsion:
        ff.MMFFAddPositionConstraint(i, 0, 1.e4)
    ff.Minimize(maxIts=10000)

    return conf
