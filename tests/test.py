import QCflow
import pytest
import sys
import unittest
import rdkit
import os
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

test_dic = {'b_18_v2' : 'COc1cc2cc(-c3nccc4nsnc34)n(C)c2cc1OC'}

#os.chdir('./tests')

def test_QCflow_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "QCflow" in sys.modules

def test_rdkit_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "rdkit" in sys.modules

def test_psi4_imported():
    assert "psi4" in sys.modules

def test_load_data():
    test_data = QCflow.load_gaussian.load_data('b_18_v2', 'opt')

    assert test_data.charge == 0

def test_open_dictionary():

    example_dic = QCflow.load_gaussian.open_dictionary('example_dic.json')

    assert example_dic['b_18_v2'] == 'COc1cc2cc(-c3nccc4nsnc34)n(C)c2cc1OC'

def test_dictionary():
    fake_dict = {'key1': 'value1', 'key2': 'value2'}
    QCflow.load_gaussian.save_dictionary(fake_dict, 'fake_dict.json')
    
    # Load the dictionary back to verify it was saved correctly
    loaded_dict = QCflow.load_gaussian.open_dictionary('fake_dict.json')
    
    assert fake_dict == loaded_dict

def test_data_dic():
    #retreave the data for all the molecules in the dictionary
    mol_data = QCflow.load_gaussian.data_dic(test_dic, 'opt')

    assert mol_data['b_18_v2'].charge == 0

def test_extract_data_from_txt():

    #gt the data from the txt file
    data = QCflow.energy_calculations.extract_data_from_txt('12/12_opt_energy_and_gap.txt')

    #check if the data is correct
    assert data['optimized_energy'] == -31659.618449
    assert data['homo'] == -5.053209
    assert data['lumo'] == -1.452872
    assert data['energy_gap'] == 3.600337

def test_cal_HOMO():
    data = QCflow.load_gaussian.load_data('b_18_v2', 'opt')

    HOMO = QCflow.energy_calculations.cal_HOMO(data)

    assert HOMO == -5.0474398129245

def test_cal_LUMO():
    data = QCflow.load_gaussian.load_data('b_18_v2', 'opt')

    LUMO = QCflow.energy_calculations.cal_LUMO(data)

    assert LUMO == -2.5766460503845

def test_cal_gap():
    data = QCflow.load_gaussian.load_data('b_18_v2', 'opt')

    GAP = QCflow.energy_calculations.cal_gap(data)

    assert GAP == 2.47079376254

def test_cal_IP_vertical():
    opt_n = QCflow.load_gaussian.load_data('b_18_v2', 'opt')
    ver_c = QCflow.load_gaussian.load_data('b_18_v2', 'sp_c')

    IP = QCflow.energy_calculations.cal_IP(opt_n, ver_c, 'vertical')

    assert IP == 6.4563109065711615

def test_cal_IP_adiabatic():
    opt_n = QCflow.load_gaussian.load_data('b_18_v2', 'opt')
    opt_c = QCflow.load_gaussian.load_data('b_18_v2', 'opt_c')

    IP = QCflow.energy_calculations.cal_IP(opt_n, opt_c, 'adiabatic')

    assert IP == 6.3088333630075795

def test_cal_EA_vertical():
    opt_n = QCflow.load_gaussian.load_data('b_18_v2', 'opt')
    ver_a = QCflow.load_gaussian.load_data('b_18_v2', 'sp_a')

    IP = QCflow.energy_calculations.cal_EA(opt_n, ver_a, 'vertical')

    assert IP == 1.0243159905294306

def test_cal_EA_adiabatic():
    opt_n = QCflow.load_gaussian.load_data('b_18_v2', 'opt')
    opt_a = QCflow.load_gaussian.load_data('b_18_v2', 'opt_a')

    IP = QCflow.energy_calculations.cal_EA(opt_n, opt_a, 'adiabatic')

    assert IP == 1.1768627430064953

def test_cal_reorg_gaussian():
    #testing if the reorganisation function works for Gaussian16
    opt_n = QCflow.load_gaussian.load_data('b_18_v2', 'opt')
    ver_c = QCflow.load_gaussian.load_data('b_18_v2', 'sp_c')
    opt_c = QCflow.load_gaussian.load_data('b_18_v2', 'opt_c')
    n_c_geo = QCflow.load_gaussian.load_data('b_18_v2', 'n_c_geo')

    reorg_val = QCflow.energy_calculations.cal_reorg(opt_n, ver_c, opt_c, n_c_geo, calculation_software='Gaussian')

    assert reorg_val == 0.29475399497459875

def test_cal_reorg_psi4():
    #testing if the reorganisation function works for psi4
    opt_n = QCflow.energy_calculations.extract_data_from_txt('12/12_opt_energy_and_gap.txt')
    ver_c = QCflow.energy_calculations.extract_data_from_txt('12/12_sp_c_energy_and_gap.txt')
    opt_c = QCflow.energy_calculations.extract_data_from_txt('12/12_opt_c_energy_and_gap.txt')
    n_c_geo = QCflow.energy_calculations.extract_data_from_txt('12/12_n_c_geo_energy_and_gap.txt')

    reorg_val = QCflow.energy_calculations.cal_reorg(opt_n, ver_c, opt_c, n_c_geo, calculation_software='Psi4')

    assert reorg_val == 0.4535870000036084

def test_adding_attch():
    mol = QCflow.fragments.adding_attach('C1=CC=CS1', find='[cH;x2]', get_rid='C([I])')

    assert (mol[0] == 'Ic1cccs1' and mol[1] == 'Ic1ccsc1') or (mol[0] == 'Ic1ccsc1' and mol[1] == 'Ic1cccs1')

def test_generate_attachment_points():
    test_mol_dic = {'test' : 'C1=CC=CS1'}

    test_mol_dic_attach = QCflow.fragments.generate_attachment_points(test_mol_dic, find='[cH;^2]', get_rid='C([I])')

    assert (test_mol_dic_attach == {'test_A': 'Ic1cccs1', 'test_B': 'Ic1ccsc1'}) or (test_mol_dic_attach == {'test_A': 'Ic1ccsc1', 'test_B': 'Ic1cccs1'})

def test_combine_structure():
    mol_1 = 'IC1=CC=CS1'
    mol_2 = 'IC1=CC=CS1'

    new_mol = QCflow.fragments.combine_structure(mol_1, mol_2)

    assert new_mol == 'c1csc(C2:cccs:2)c1'

def test_make_molecule_dic_from_2_dic():
    test_dic_1 = {'mol_A' : 'IC1=CC=CS1'}
    test_dic_2 = {'mol_B' : 'IC1=CC=CS1'}

    new_dic = QCflow.fragments.make_molecule_dic_from_2_dic(test_dic_1, test_dic_2)

    assert new_dic == {'mol_A_mol_B' : 'c1csc(C2:cccs:2)c1'}

def test_getBond():
    tor_mol = Chem.MolFromSmiles('COc1cc2cc(-c3nccc4nsnc34)n(C)c2cc1OC')

    bond = QCflow.find_torsion.getBond(tor_mol)

    assert bond[0] == (6, 7)

def test_getTorsion():
    tor_mol = Chem.MolFromSmiles('COc1cc2cc(-c3nccc4nsnc34)n(C)c2cc1OC')

    bond = QCflow.find_torsion.getBond(tor_mol)

    torsion = QCflow.find_torsion.getTorsion(tor_mol, bond[0])

    assert torsion == (16, 6, 7, 8)

def test_embed_molecule():
    tor_mol = Chem.MolFromSmiles('COc1cc2cc(-c3nccc4nsnc34)n(C)c2cc1OC')

    mol3d = QCflow.find_torsion.embed_molecule(tor_mol)

    assert Chem.MolToSmiles(mol3d) == '[H]c1nc(-c2c([H])c3c([H])c(OC([H])([H])[H])c(OC([H])([H])[H])c([H])c3n2C([H])([H])[H])c2nsnc2c1[H]'

def test_success_test():
    test_dic = {'b_18_v2' : 'COc1cc2cc(-c3nccc4nsnc34)n(C)c2cc1OC'}

    test_dic_data = {}
    for k, v in test_dic.items():
        test_dic_data[k] = QCflow.load_gaussian.load_data(k, 'opt')

    passed, failed = QCflow.testing_data.success_test(test_dic_data)

    assert len(passed) == 1 and len(failed) == 0

def test_functional_test():
    test_dic = {'b_18_v2' : 'COc1cc2cc(-c3nccc4nsnc34)n(C)c2cc1OC'}

    test_dic_data = {}
    for k, v in test_dic.items():
        test_dic_data[k] = QCflow.load_gaussian.load_data(k, 'opt')

    passed, failed = QCflow.testing_data.functional_test(test_dic_data, functional='B3LYP')

    assert len(passed) == 1 and len(failed) == 0

def test_basis_set_test():
    test_dic = {'b_18_v2' : 'COc1cc2cc(-c3nccc4nsnc34)n(C)c2cc1OC'}

    test_dic_data = {}
    for k, v in test_dic.items():
        test_dic_data[k] = QCflow.load_gaussian.load_data(k, 'opt')

    passed, failed = QCflow.testing_data.basis_set_test(test_dic_data, basis_set='6-31G(d)')

    assert len(passed) == 1 and len(failed) == 0

def test_find_min_energy_index():
    tor_1_data = QCflow.load_gaussian.load_data('torsion_1', 'tor')

    tor_index = QCflow.torsion_parser.find_min_energy_index(tor_1_data)

    assert tor_index == 18

def test_min_angle():
    tor_1_data = QCflow.load_gaussian.load_data('torsion_1', 'tor')

    min_ang = QCflow.torsion_parser.min_angle(tor_1_data)

    assert min_ang == -179.9988

def test_torsion_parser():
    torsion_smi = 'COC(C(OC)=C1)=CC2=C1C=C(N2[H])C3=C(OC)C=CS3'
    
    result = QCflow.torsion_parser.torsion_parser('torsion_1', torsion_smi)

    assert isinstance(result, Chem.rdchem.Conformer)  # Assert that result is an rdkit.Chem.rdchem.Conformer

def test_setting_dihedral():
    torsion_smi = 'COC(C(OC)=C1)=CC2=C1C=C(N2[H])C3=C(OC)C=CS3'
    
    result = QCflow.torsion_parser.setting_dihedral(torsion_smi, 90)

    mol = Chem.MolFromSmiles(torsion_smi)
    #Converts the SMILE into rdkit readable string

    bond = QCflow.find_torsion.getBond(mol)

    torsion_atoms = QCflow.find_torsion.getTorsion(mol, bond[0])

    angle = rdMolTransforms.GetDihedralDeg(result, torsion_atoms[0], torsion_atoms[1], torsion_atoms[2], torsion_atoms[3])

    assert isinstance(result, Chem.rdchem.Conformer)  and int(angle) == 90

def test_find_dihedral_opt():
    opt_data = QCflow.load_gaussian.load_data('b_18_v2', 'opt')

    opt_smi = 'COc1cc2cc(-c3nccc4nsnc34)n(C)c2cc1OC'

    result = QCflow.torsion_parser.finding_dihedral_opt(opt_smi, opt_data)

    assert result == 9.918549403987761

def test_find_planarity():
    angle = 90
    result = QCflow.torsion_parser.find_planarity(angle)
    assert int(result) == 0

def test_getBondLinkers():
    linker_type = {
    'mol_1' : 'single',
    'mol_2' : 'double',
    'mol_3' : 'imine',
    'mol_4' : 'thio'
    }

    mol_smi = {
    'mol_1' : 'COC1=CC2=C(C=C1OC)N(C)C(C3=NC=CC4=NSN=C34)=C2',
    'mol_2' : 'COC1=CC2=C(C=C1OC)N(C)C(/C=C/C3=NC=CC4=NSN=C34)=C2',
    'mol_3' : 'COC1=CC2=C(C=C1OC)N(C)C(/N=C/C3=NC=CC4=NSN=C34)=C2',
    'mol_4' : 'COC1=CC2=C(C=C1OC)N(C)C(C(S3)=CC=C3C4=NC=CC5=NSN=C45)=C2'
    }

    bonds_link = {}
    for k, v in linker_type.items():
        mol = Chem.MolFromSmiles(mol_smi[k])
        bond = QCflow.torsion_parser.getBondLinkers(mol, v)
        bonds_link[k] = bond

    assert bonds_link == {'mol_1': ((12, 13),), 'mol_2': ((12, 13, 14, 15),), 'mol_3': ((12, 13, 14, 15),), 'mol_4': ((12, 13), (17, 18))}

def test_getTorsion_one():
    mol_2_mol = Chem.MolFromSmiles('COC1=CC2=C(C=C1OC)N(C)C(/C=C/C3=NC=CC4=NSN=C34)=C2')

    bonds = ((12, 13, 14, 15),)

    result = QCflow.torsion_parser.getTorsion_one(mol_2_mol, bonds[0])

    assert result == (10, 12, 13, 14)

def test_getTorsion_two():
    mol_2_mol = Chem.MolFromSmiles('COC1=CC2=C(C=C1OC)N(C)C(/C=C/C3=NC=CC4=NSN=C34)=C2')

    bonds = ((12, 13, 14, 15),)

    result = QCflow.torsion_parser.getTorsion_two(mol_2_mol, bonds[0])

    assert result == (13, 14, 15, 16)

def test_write_gaussian():
    mol_smi = 'COc1cc2cc(-c3nccc4nsnc34)n(C)c2cc1OC'
    mol = Chem.MolFromSmiles(mol_smi)
    mol3d = QCflow.find_torsion.embed_molecule(mol)
    conf = mol3d.GetConformer()

    os.chdir('b_18_v2')

    QCflow.write_gaussian.write_gaussian('opt', 'b_18_v2', mol_smi, functional='B3LYP', basis_set='6-31G*', mol=mol3d, torsion=0, conformer=conf)

    # Check if the .com file exists
    file_path = 'b_18_v2_opt.com'
    assert os.path.exists(file_path), f"{file_path} does not exist."

    # Check if the .com file contains the correct information
    with open(file_path, 'r') as file:
        content = file.read()
        assert 'B3LYP' in content, "Functional B3LYP not found in the file."
        assert '6-31G*' in content, "Basis set 6-31G* not found in the file."
        assert '0 1' in content, "Molecule charge and multiplicity not found in the file."
        assert 'b_18_v2 opt Smile String: COc1cc2cc(-c3nccc4nsnc34)n(C)c2cc1OC' in content, "Molecule name and SMILES string not found in the file."

    os.chdir('../')

def test_write_slurm():
    os.chdir('b_18_v2')

    QCflow.slurm.write_slurm('opt', 'b_18_v2', 20)

    # Check if the .slurm file exists
    file_path = 'b_18_v2_opt.sh'
    assert os.path.exists(file_path), f"{file_path} does not exist."

    # Check if the .slurm file contains the correct information
    with open(file_path, 'r') as file:
        content = file.read()
        assert '#SBATCH --job-name=b_18_v2_opt' in content, "Job name not found in the file."
        assert '#SBATCH --ntasks=20' in content, "Number of CPUs  not found in the file."
        assert '#SBATCH -p cpu ' in content, "cpu partition  not found in the file."
        assert 'g16 $INPUTFILE > $OUTPUTFILE' in content, "g16 execution command not found in the file."

    os.chdir('../')

def test_write_slurm_psi4():
    os.chdir('1')

    QCflow.slurm.write_slurm_psi4('opt', '1', time=24, cpus=10)

    # Check if the .slurm file exists
    file_path = '1_opt.sh'
    assert os.path.exists(file_path), f"{file_path} does not exist."

    # Check if the .slurm file contains the correct information
    with open(file_path, 'r') as file:
        content = file.read()
        assert '#SBATCH --job-name=1_opt' in content, "Job name not found in the file."
        assert '#SBATCH --ntasks=10' in content, "Number of CPUs  not found in the file."
        assert '#SBATCH -p cpu ' in content, "cpu partition  not found in the file."
        assert 'module load cuda/10.0.130-gcc-13.2.0' in content, "cuda module load command not found in the file."
    
    os.chdir('../')

def test_rdkit_predict_conf():
    test_smi = 'COc1cc2cc(-c3nccc4nsnc34)n(C)c2cc1OC'
    conf = QCflow.run_gaussian.rdkit_predict_conf(test_smi)

    assert isinstance(conf, Chem.rdchem.Conformer)

def test_write_psi4():
    smi = 'CC(=O)C1=CC=C(C=C1)C(=O)C'
    mol_name = '1'
    mol = Chem.MolFromSmiles(smi)
    #gets rdkit estimated coordinates of dimer
    mol3d = QCflow.find_torsion.embed_molecule(mol)
    conf_geo = QCflow.run_gaussian.rdkit_predict_conf(smi)
    #go into the directory
    os.chdir('1')
    #write the input file
    QCflow.write_psi4.write_psi4('opt', mol_name, smi, functional='b3lyp', basis_set='6-31g*', mol=mol3d, conformer=conf_geo)
    #go back to the original directory
    os.chdir('../')
    file_path = '1/1_opt.py'
    assert os.path.exists(file_path), 'The file was not created'
    with open(file_path, 'r') as file:
        content = file.read()
        assert 'import psi4' in content, "Psi4 import not written in the file"
        assert 'psi4.set_options' in content, "Psi4 options not written in the file"
        assert 'psi4.set_memory' in content, "Psi4 memory not written in the file"
        assert 'optimized_geometry_xyz' in content, "Optimized geometry not written in the file"

def test_write_psi4_reorg_test():
    mol_name = '1'
    #go into the directory
    os.chdir('1')
    #write the input file
    QCflow.write_psi4.write_psi4_reorg('cation', mol_name, functional='b3lyp', basis_set='6-31g*')
    #go back to the original directory
    os.chdir('../')
    file_path = '1/1_cation.py'
    assert os.path.exists(file_path), 'The file was not created'
    with open(file_path, 'r') as file:
        content = file.read()
        assert 'import psi4' in content, "Psi4 import not written in the file"
        assert 'psi4.set_options' in content, "Psi4 options not written in the file"
        assert 'psi4.set_memory' in content, "Psi4 memory not written in the file"
        assert 'optimized_geometry_xyz' in content, "Optimized geometry not written in the file"
