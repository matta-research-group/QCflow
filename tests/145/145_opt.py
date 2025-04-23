import psi4
 
psi4.set_memory('4000 mb')
psi4.set_num_threads(4)
 
mol_name = psi4.geometry("""
 0 1 
 C -7.75892 -0.52073 1.91815 
 O -6.66748 0.32284 2.27720 
 C -5.64012 0.32547 1.36289 
 C -4.36788 0.16124 1.93207 
 C -3.23338 0.15693 1.10327 
 C -1.83762 0.00380 1.33314 
 C -1.18365 0.08473 0.11708 
 C 0.21688 -0.00540 -0.22003 
 S 1.39883 -0.26306 1.00087 
 C 2.68188 -0.26212 -0.12469 
 C 4.05967 -0.46630 0.34866 
 C 4.64191 -1.75539 0.35611 
 C 3.85322 -2.90059 -0.15193 
 O 4.38612 -3.83641 -0.73376 
 C 5.95791 -1.91753 0.80262 
 C 6.63168 -3.23371 0.92046 
 O 7.81973 -3.30845 1.21907 
 C 6.68428 -0.79301 1.21379 
 C 6.13637 0.47749 1.19624 
 N 6.90095 1.50669 1.60377 
 C 6.34085 2.73310 1.57486 
 C 5.04990 2.93072 1.14553 
 N 4.27223 1.90818 0.73088 
 C 4.80266 0.66062 0.75742 
 C 2.23118 -0.05637 -1.41154 
 C 0.80999 0.09171 -1.47250 
 O 0.03841 0.29992 -2.57704 
 C 0.75814 0.44077 -3.79561 
 N -2.15673 0.28724 -0.83475 
 C -3.40713 0.33533 -0.26470 
 C -4.66587 0.53075 -0.85353 
 C -5.79457 0.54096 -0.01769 
 O -7.07036 0.76461 -0.47206 
 C -7.22667 1.18292 -1.82103 
 H -8.53166 0.05830 1.40432 
 H -7.45072 -1.38311 1.31582 
 H -8.19837 -0.90468 2.84418 
 H -4.26723 0.02779 3.00490 
 H -1.37150 -0.14986 2.29761 
 H 2.75964 -2.87239 -0.02150 
 H 6.04469 -4.15045 0.77006 
 H 7.71030 -0.90535 1.56614 
 H 6.96566 3.55354 1.90838 
 H 4.59869 3.91592 1.12138 
 H 2.90219 -0.01192 -2.25995 
 H 0.03164 0.62309 -4.59335 
 H 1.30017 -0.47829 -4.04243 
 H 1.43555 1.30062 -3.75871 
 H -1.96102 0.38985 -1.82166 
 H -4.73473 0.67344 -1.92448 
 H -8.28654 1.40321 -1.98120 
 H -6.94626 0.38342 -2.51433 
 H -6.66294 2.09992 -2.02248 
""") 
 
psi4.set_module_options('scf', {'maxiter': 75}) 
psi4.set_options({'basis': '6-31g*', 'scf_type': 'df'})
 
energy, wfn = psi4.optimize('b3lyp', engine='geometric', return_wfn=True)  
 
optimized_geometry_xyz = wfn.molecule().to_string('xyz') 
with open('145_opt.xyz', 'w') as xyz_file:
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
 
with open('145_opt_energy_and_gap.txt', 'w') as file:
    file.write(f"Optimized energy: {energy_ev:.6f} eV\n") 
    file.write(f"HOMO: {homo_energy_ev:.6f} eV\n") 
    file.write(f"LUMO: {lumo_energy_ev:.6f} eV\n") 
    file.write(f"Energy gap (HOMO-LUMO): {energy_gap_ev:.6f} eV\n") 
 
