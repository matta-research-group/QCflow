import psi4
 
psi4.set_memory('4000 mb')
psi4.set_num_threads(4)
 
mol = psi4.geometry(open('1_opt.xyz').read())
psi4.core.set_active_molecule(mol)
mol.set_molecular_charge(1) 
mol.set_multiplicity(2) 
 
psi4.set_module_options('scf', {'maxiter': 100}) 
psi4.set_options({'basis': '6-31g*', 'scf_type': 'df', 'reference': 'uhf'})
 
energy, wfn = psi4.optimize('b3lyp', return_wfn=True)  
 
optimized_geometry_xyz = wfn.molecule().to_string('xyz') 
with open('1_opt_c.xyz', 'w') as xyz_file:
    xyz_file.write(optimized_geometry_xyz) 
 
psi4.core.set_active_molecule(wfn.molecule()) 
 
orbital_energies = wfn.epsilon_a_subset("AO", "ALL").np 
 
n_occ_alpha = wfn.nalpha() 
homo_energy = orbital_energies[n_occ_alpha - 1] 
lumo_energy = orbital_energies[n_occ_alpha] 
energy_gap = lumo_energy - homo_energy 
 
hartree_to_ev = 27.2114079527 
energy_ev =  energy * hartree_to_ev
homo_energy_ev = homo_energy * hartree_to_ev 
lumo_energy_ev = lumo_energy * hartree_to_ev 
energy_gap_ev = lumo_energy_ev -  homo_energy_ev 
 
with open('1_opt_c_energy_and_gap.txt', 'w') as file:
    file.write(f"Optimized energy: {energy_ev:.6f} eV\n") 
    file.write(f"HOMO: {homo_energy_ev:.6f} eV\n") 
    file.write(f"LUMO: {lumo_energy_ev:.6f} eV\n") 
    file.write(f"Energy gap (HOMO-LUMO): {energy_gap_ev:.6f} eV\n") 
 
 
mol_sp = psi4.geometry(open('1_opt_c.xyz').read())
psi4.core.set_active_molecule(mol_sp)
mol_sp.set_molecular_charge(0) 
mol_sp.set_multiplicity(1) 
 
psi4.set_options({'basis': '6-31g*', 'scf_type': 'df'})
 
energy, wfn = psi4.energy('b3lyp', return_wfn=True)  
 
 
psi4.core.set_active_molecule(wfn.molecule()) 
 
orbital_energies = wfn.epsilon_a_subset("AO", "ALL").np 
 
n_occ_alpha = wfn.nalpha() 
homo_energy = orbital_energies[n_occ_alpha - 1] 
lumo_energy = orbital_energies[n_occ_alpha] 
energy_gap = lumo_energy - homo_energy 
 
energy_ev =  energy * hartree_to_ev
homo_energy_ev = homo_energy * hartree_to_ev 
lumo_energy_ev = lumo_energy * hartree_to_ev 
energy_gap_ev = lumo_energy_ev -  homo_energy_ev 
 
with open('1_n_c_geo_energy_and_gap.txt', 'w') as file:
    file.write(f"Single Point energy: {energy_ev:.6f} eV\n") 
    file.write(f"HOMO: {homo_energy_ev:.6f} eV\n") 
    file.write(f"LUMO: {lumo_energy_ev:.6f} eV\n") 
    file.write(f"Energy gap (HOMO-LUMO): {energy_gap_ev:.6f} eV\n") 
 
