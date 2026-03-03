from QCflow.fragments import *
from QCflow.find_torsion import *
from QCflow.write_gaussian import *
from QCflow.slurm import *
from QCflow.torsion_parser import *
from QCflow.run_gaussian import *
from QCflow.write_psi4 import *
import os
import json
import subprocess
import tempfile

def run_xtb_optimisation(xyz: str, output_file: str, charge=0, uhf=0, gfn=2):
    """
    Runs a xTB geometry optimisation using GFN2-xTB and saves the result to output_file.
    
    Parameters
    ----------
    xyz (str): The input geometry in XYZ format as a string.
    output_file (str): The path to the file where the optimised geometry will be saved
    charge (int, optional): The total charge of the system. Default is 0.
    uhf (int, optional): The number of unpaired electrons (spin multiplicity - 1). Default is 0 (singlet).
    gfn (int, optional): The GFN-xTB method to use (1, 2, or 3). Default is 2 (GFN2-xTB).
    
    Returns
    -------
         None: Writes the optimised geometry to output_file.

    """
    with tempfile.TemporaryDirectory() as tmpdir:
        xyz_path = os.path.join(tmpdir, "input.xyz")
        with open(xyz_path, 'w') as f:
            f.write(xyz)

        cmd = [
            "xtb", xyz_path,
            "--opt",
            "--chrg", str(charge),
            "--uhf", str(uhf),
            "--gfn", str(gfn)
        ]

        try:
            subprocess.run(cmd, cwd=tmpdir, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"xTB failed: {e.stderr.decode()}")

        opt_path = os.path.join(tmpdir, "xtbopt.xyz")
        if not os.path.exists(opt_path):
            raise FileNotFoundError("Optimised structure (xtbopt.xyz) not found.")

        with open(opt_path, 'r') as f:
            optimised_xyz = f.read()

        with open(output_file, 'w') as f:
            f.write(optimised_xyz)

        print(f"Optimised structure saved to: {output_file}")

def load_xyz_from_file(filepath):
    """
    Loads an XYZ file and returns the content as a string.
    
    Parameters
    ----------
    filepath (str): The path to the XYZ file to be loaded.
    
    Returns
    -------
         str: The content of the XYZ file as a string.
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Remove the first two lines (atom count and comment/energy)
    coordinates = lines[2:]

    # Recalculate number of atoms just in case
    num_atoms = len(coordinates)

    # Build the xyz string
    xyz_string = f"{num_atoms}\n\n" + "".join(coordinates)
    return xyz_string


def write_xTB_psi4(job_name, mol_name, smile, functional='b3lyp', basis_set='6-31g*', mol=None, conformer=None):
    """
    Runs a xTB geometry optimisation and caulcates electronic proeprities using Psi4, based on the provided parameters.
    
    Parameters
    ----------
    job_name (str): The type of job to run. Possible values:
        - 'opt': Optimisation neutral                                                                                                                               
    mol_name (str): The name of the molecule.
    smile (str): The SMILE string of the molecule.
    functional (str, optional): The functional to use. Default is 'b3lyp'.
    basis_set (str, optional): The basis set to use. Default is '6-31g'.
    mol (rdkit.Chem.rdchem.Mol, optional): The RDKit embedded molecule object.
    conformer (rdkit.Chem.rdchem.Conformer, optional): The RDKit conformer of the molecule.
    
    Returns
    -------
        None: Writes the psi4 input file to disk.

    """
    if (job_name=='opt'):

        old_chk = ' '
        symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        # get x y z coords
        geometry = conformer
        calculation = 'opt'
        torsion_data = f' \n'

        if Descriptors.NumRadicalElectrons(mol) == 0:
            mult_chg = '0 1'
        if Descriptors.NumRadicalElectrons(mol) == 1:
            mult_chg = '0 2'
        if Descriptors.NumRadicalElectrons(mol) == 2:
            mult_chg = '0 3' 
                                                                                          
    file_name = f'{mol_name}_{job_name}_xTB.py'
    #Includes information about basis set to allow for ramping
    title = f'{mol_name} {job_name} Smile String: {smile}'
    number_proc = 'set_num_threads(4)'
    memory_num = '4000 mb'
    number_of_atoms = mol.GetNumAtoms()

    with open(file_name, 'w') as file:
        file.write('import subprocess\n')
        file.write('import tempfile\n')
        file.write('import os\n')
        file.write('import psi4\n')
        file.write('from QCflow.write_xTB import *\n')
        file.write(' \n')#
        file.write(f'xyz_string = """{number_of_atoms}\n')
        file.write(' \n')
        if old_chk==' ':
            for atom,symbol in enumerate(symbols):
                p = geometry.GetAtomPosition(atom)
                # atom  x y z
                line = f' {symbol} {p.x:.5f} {p.y:.5f} {p.z:.5f} \n'
                file.write(line)
        file.write('""" \n')
        file.write(' \n')
        file.write('run_xtb_optimisation(xyz_string, "optimised_structure.xyz", gfn=2)\n')
        file.write(' \n')
        file.write(f"mol_sp = psi4.geometry(open('optimised_structure.xyz').read())\n")#use geometry from opt
        file.write('psi4.core.set_active_molecule(mol_sp)\n')#set as active molecule
        file.write(f"mol_sp.set_molecular_charge(0) \n")# set charge
        file.write(f"mol_sp.set_multiplicity(1) \n")# set multiplicity
        file.write(' \n')
        file.write(f"psi4.set_options({{'basis': '{basis_set}', 'scf_type': 'df'}})\n")
        file.write(' \n')
        file.write(f"energy, wfn = psi4.energy('{functional}', return_wfn=True)  \n")
        file.write(' \n')
        file.write(' \n')
        file.write(f"psi4.core.set_active_molecule(wfn.molecule()) \n")
        file.write(' \n')
        file.write('orbital_energies = wfn.epsilon_a_subset("AO", "ALL").np \n')
        file.write(' \n')
        file.write('n_occ_alpha = wfn.nalpha() \n')
        file.write('homo_energy = orbital_energies[n_occ_alpha - 1] \n')
        file.write('lumo_energy = orbital_energies[n_occ_alpha] \n')
        file.write('energy_gap = lumo_energy - homo_energy \n')
        file.write(' \n')
        file.write('hartree_to_ev = 27.2114079527\n')
        file.write('energy_ev =  energy * hartree_to_ev\n')
        file.write('homo_energy_ev = homo_energy * hartree_to_ev \n')
        file.write('lumo_energy_ev = lumo_energy * hartree_to_ev \n')
        file.write('energy_gap_ev = lumo_energy_ev -  homo_energy_ev \n')
        file.write(' \n')
        file.write(f"with open('{mol_name}_xTB_energy_and_gap.txt', 'w') as file:\n")
        file.write('    file.write(f"Single Point energy: {energy_ev:.6f} eV\\n") \n')
        file.write('    file.write(f"HOMO: {homo_energy_ev:.6f} eV\\n") \n')
        file.write('    file.write(f"LUMO: {lumo_energy_ev:.6f} eV\\n") \n')
        file.write('    file.write(f"Energy gap (HOMO-LUMO): {energy_gap_ev:.6f} eV\\n") \n')
        file.write(' \n')


def write_xTB_psi4_reorg(job_name, mol_name, functional='b3lyp', basis_set='6-31g*'):
    """
    Generates a xTB and psi4 input file based on the provided parameters for reorganisation calculations.
    
    Parameters
    ----------
    job_name (str): The type of job to run. Possible values:
        - 'cation': Geometry optimisation cation (opt_c) and single of neutral charge, cation geometry (n_c_geo)
        - 'anion': Geometry optimisation anion (opt_a) and single of neutral charge, anion geometry (n_a_geo)
        - 'sp_c': Single point calculation of neutral geometry at cation charge
        - 'sp_a': Single point calculation of neutral geometry at anion charge                                                                                                                                     
    mol_name (str): The name of the molecule.
    functional (str, optional): The functional to use. Default is 'b3lyp'.
    basis_set (str, optional): The basis set to use. Default is '6-31g'.
    
    Returns
    -------
         None: Writes the psi4 input file to disk.

    """
    if (job_name=='cation'):

        charge = 1
        mult = 2

        job_type = 'opt_c'

        job_sp = 'n_c_geo'

        with open('optimised_structure.xyz', 'r') as f:
            lines = f.readlines()
            number_of_atoms = int(lines[0].strip())
            coordinate_lines = [line.strip() for line in lines[2:2 + number_of_atoms]]


    if (job_name=='anion'):

        charge = -1
        mult = 2

        job_type = 'opt_a'

        job_sp = 'n_a_geo'

        with open('optimised_structure.xyz', 'r') as f:
            lines = f.readlines()
            number_of_atoms = int(lines[0].strip())
            coordinate_lines = [line.strip() for line in lines[2:2 + number_of_atoms]]

    if (job_name=='sp_c'):

        charge = 1
        mult = 2

        with open('optimised_structure.xyz', 'r') as f:
            lines = f.readlines()
            number_of_atoms = int(lines[0].strip())
            coordinate_lines = [line.strip() for line in lines[2:2 + number_of_atoms]]

    if (job_name=='sp_a'):

        charge = -1
        mult = 2

        with open('optimised_structure.xyz', 'r') as f:
            lines = f.readlines()
            number_of_atoms = int(lines[0].strip())
            coordinate_lines = [line.strip() for line in lines[2:2 + number_of_atoms]]

                                                                                          
    file_name = f'{mol_name}_{job_name}.py'
    #Includes information about basis set to allow for ramping
    memory_num = '4000 mb'

    if (job_name=='anion') or (job_name=='cation'):

        with open(file_name, 'w') as file:
            file.write('import subprocess\n')
            file.write('import tempfile\n')
            file.write('import os\n')
            file.write('import psi4\n')
            file.write('from QCflow.write_xTB import *\n')
            file.write(' \n')#
            file.write(f'xyz_string = """{number_of_atoms}\n')
            file.write(' \n')
            # Now insert the parsed coordinate lines
            for line in coordinate_lines:
                file.write(f'{line}\n')
            file.write('""" \n')
            file.write(' \n')
            file.write(f"run_xtb_optimisation(xyz_string, '{job_type}_optimised_structure.xyz', charge={charge}, uhf={mult}, gfn=2)\n")
            file.write(' \n')
            file.write(f"mol = psi4.geometry(open('{job_type}_optimised_structure.xyz').read())\n")#use geometry from opt
            file.write('psi4.core.set_active_molecule(mol)\n')#set as active molecule
            file.write(f"mol.set_molecular_charge({charge}) \n")# set charge
            file.write(f"mol.set_multiplicity({mult}) \n")# set multiplicity
            file.write("psi4.set_module_options('scf', {'maxiter': 100}) \n")
            file.write(f"psi4.set_options({{'basis': '{basis_set}', 'scf_type': 'df', 'reference': 'uhf'}})\n")
            file.write(' \n')
            file.write(f"energy, wfn = psi4.optimize('{functional}', engine='geometric', return_wfn=True)  \n")
            file.write(' \n')
            file.write("optimized_geometry_xyz = wfn.molecule().to_string('xyz') \n")
            file.write(f"with open('{mol_name}_{job_type}.xyz', 'w') as xyz_file:\n")
            file.write("    xyz_file.write(optimized_geometry_xyz) \n")
            file.write(' \n')
            file.write(f"psi4.core.set_active_molecule(wfn.molecule()) \n")
            file.write(' \n')
            file.write('orbital_energies = wfn.epsilon_a_subset("AO", "ALL").np \n')
            file.write(' \n')
            file.write('n_occ_alpha = wfn.nalpha() \n')
            file.write('homo_energy = orbital_energies[n_occ_alpha - 1] \n')
            file.write('lumo_energy = orbital_energies[n_occ_alpha] \n')
            file.write('energy_gap = lumo_energy - homo_energy \n')
            file.write(' \n')
            file.write('hartree_to_ev = 27.2114079527 \n')
            file.write('energy_ev =  energy * hartree_to_ev\n')
            file.write('homo_energy_ev = homo_energy * hartree_to_ev \n')
            file.write('lumo_energy_ev = lumo_energy * hartree_to_ev \n')
            file.write('energy_gap_ev = lumo_energy_ev -  homo_energy_ev \n')
            file.write(' \n')
            file.write(f"with open('{mol_name}_{job_type}_energy_and_gap.txt', 'w') as file:\n")
            file.write('    file.write(f"Optimized energy: {energy_ev:.6f} eV\\n") \n')
            file.write('    file.write(f"HOMO: {homo_energy_ev:.6f} eV\\n") \n')
            file.write('    file.write(f"LUMO: {lumo_energy_ev:.6f} eV\\n") \n')
            file.write('    file.write(f"Energy gap (HOMO-LUMO): {energy_gap_ev:.6f} eV\\n") \n')
            file.write(' \n')# Now do the n_c_geo or n_a_geo part
            file.write(f"mol_sp = psi4.geometry(open('{mol_name}_{job_type}.xyz').read())\n")#use geometry from opt
            file.write('psi4.core.set_active_molecule(mol_sp)\n')#set as active molecule
            file.write(f"mol_sp.set_molecular_charge(0) \n")# set charge
            file.write(f"mol_sp.set_multiplicity(1) \n")# set multiplicity
            file.write(' \n')
            file.write(f"psi4.set_options({{'basis': '{basis_set}', 'scf_type': 'df'}})\n")
            file.write(' \n')
            file.write(f"energy, wfn = psi4.energy('{functional}', return_wfn=True)  \n")
            file.write(' \n')
            file.write(' \n')
            file.write(f"psi4.core.set_active_molecule(wfn.molecule()) \n")
            file.write(' \n')
            file.write('orbital_energies = wfn.epsilon_a_subset("AO", "ALL").np \n')
            file.write(' \n')
            file.write('n_occ_alpha = wfn.nalpha() \n')
            file.write('homo_energy = orbital_energies[n_occ_alpha - 1] \n')
            file.write('lumo_energy = orbital_energies[n_occ_alpha] \n')
            file.write('energy_gap = lumo_energy - homo_energy \n')
            file.write(' \n')
            file.write('energy_ev =  energy * hartree_to_ev\n')
            file.write('homo_energy_ev = homo_energy * hartree_to_ev \n')
            file.write('lumo_energy_ev = lumo_energy * hartree_to_ev \n')
            file.write('energy_gap_ev = lumo_energy_ev -  homo_energy_ev \n')
            file.write(' \n')
            file.write(f"with open('{mol_name}_{job_sp}_energy_and_gap.txt', 'w') as file:\n")
            file.write('    file.write(f"Single Point energy: {energy_ev:.6f} eV\\n") \n')
            file.write('    file.write(f"HOMO: {homo_energy_ev:.6f} eV\\n") \n')
            file.write('    file.write(f"LUMO: {lumo_energy_ev:.6f} eV\\n") \n')
            file.write('    file.write(f"Energy gap (HOMO-LUMO): {energy_gap_ev:.6f} eV\\n") \n')
            file.write(' \n')

        if (job_name=='sp_c') or (job_name=='sp_a'):
            with open(file_name, 'w') as file:
                file.write('import psi4\n')
                file.write(' \n')#
                file.write(f"psi4.set_memory('{memory_num}')\n")
                file.write('psi4.set_num_threads(4)\n')
                file.write(' \n')#
                file.write(f"mol = psi4.geometry(open('optimised_structure.xyz').read())\n")#use geometry from opt
                file.write('psi4.core.set_active_molecule(mol)\n')#set as active molecule
                file.write(f"mol.set_molecular_charge({charge}) \n")# set charge
                file.write(f"mol.set_multiplicity({mult}) \n")# set multiplicity
                file.write(' \n')
                file.write("psi4.set_module_options('scf', {'maxiter': 50}) \n")
                file.write(f"psi4.set_options({{'basis': '{basis_set}', 'scf_type': 'df', 'reference': 'uhf'}})\n")
                file.write(' \n')
                file.write(f"energy, wfn = psi4.energy('{functional}', return_wfn=True)  \n")
                file.write(' \n')
                file.write(' \n')
                file.write(f"psi4.core.set_active_molecule(wfn.molecule()) \n")
                file.write(' \n')
                file.write('orbital_energies = wfn.epsilon_a_subset("AO", "ALL").np \n')
                file.write(' \n')
                file.write('n_occ_alpha = wfn.nalpha() \n')
                file.write('homo_energy = orbital_energies[n_occ_alpha - 1] \n')
                file.write('lumo_energy = orbital_energies[n_occ_alpha] \n')
                file.write('energy_gap = lumo_energy - homo_energy \n')
                file.write(' \n')
                file.write('hartree_to_ev = 27.2114079527 \n')
                file.write('energy_ev =  energy * hartree_to_ev\n')
                file.write('homo_energy_ev = homo_energy * hartree_to_ev \n')
                file.write('lumo_energy_ev = lumo_energy * hartree_to_ev \n')
                file.write('energy_gap_ev = lumo_energy_ev -  homo_energy_ev \n')
                file.write(' \n')
                file.write(f"with open('{mol_name}_{job_name}_energy_and_gap.txt', 'w') as file:\n")
                file.write('    file.write(f"Single Point energy: {energy_ev:.6f} eV\\n") \n')
                file.write('    file.write(f"HOMO: {homo_energy_ev:.6f} eV\\n") \n')
                file.write('    file.write(f"LUMO: {lumo_energy_ev:.6f} eV\\n") \n')
                file.write('    file.write(f"Energy gap (HOMO-LUMO): {energy_gap_ev:.6f} eV\\n") \n')
                file.write(' \n')
