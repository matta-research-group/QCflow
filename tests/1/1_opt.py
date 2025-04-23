import psi4
 
psi4.set_memory('4000 mb')
psi4.set_num_threads(4)
 
mol_name = psi4.geometry("""
 0 1 
 C 3.85817 0.25578 -0.27077 
 C 2.64161 0.23270 -1.15846 
 O 2.74080 0.34537 -2.37861 
 C 1.29339 0.06954 -0.53697 
 C 0.16945 0.04976 -1.37512 
 C -1.11711 -0.10063 -0.84395 
 C -1.29339 -0.23335 0.53697 
 C -0.16945 -0.21357 1.37512 
 C 1.11711 -0.06318 0.84395 
 C -2.64161 -0.39651 1.15846 
 O -2.74080 -0.50917 2.37861 
 C -3.85817 -0.41958 0.27077 
 H 3.95113 -0.68989 0.26941 
 H 4.75205 0.38168 -0.88980 
 H 3.80533 1.09898 0.42272 
 H 0.29085 0.15205 -2.45319 
 H -1.95755 -0.11016 -1.53212 
 H -0.29085 -0.31586 2.45319 
 H 1.95755 -0.05365 1.53212 
 H -3.80533 -1.26279 -0.42272 
 H -4.75205 -0.54549 0.88980 
 H -3.95113 0.52608 -0.26940 
""") 
 
psi4.set_module_options('scf', {'maxiter': 75}) 
psi4.set_options({'basis': '6-31g*', 'scf_type': 'df'})
 
energy, wfn = psi4.optimize('b3lyp', return_wfn=True)  
 
optimized_geometry_xyz = wfn.molecule().to_string('xyz') 
with open('1_opt.xyz', 'w') as xyz_file:
    xyz_file.write(optimized_geometry_xyz) 
 
psi4.core.set_active_molecule(wfn.molecule()) 
 
orbital_energies = wfn.epsilon_a_subset("AO", "ALL").np 
 
n_occ_alpha = wfn.nalpha() 
homo_energy = orbital_energies[n_occ_alpha - 1] 
lumo_energy = orbital_energies[n_occ_alpha] 
energy_gap = lumo_energy - homo_energy 
 
hartree_to_ev = 27.2114079527 
energy_ev = energy * hartree_to_ev 
homo_energy_ev = homo_energy * hartree_to_ev 
lumo_energy_ev = lumo_energy * hartree_to_ev 
energy_gap_ev = lumo_energy_ev -  homo_energy_ev 
 
with open('1_opt_energy_and_gap.txt', 'w') as file:
    file.write(f"Optimized energy: {energy_ev:.6f} eV\n") 
    file.write(f"HOMO: {homo_energy_ev:.6f} eV\n") 
    file.write(f"LUMO: {lumo_energy_ev:.6f} eV\n") 
    file.write(f"Energy gap (HOMO-LUMO): {energy_gap_ev:.6f} eV\n") 
 
