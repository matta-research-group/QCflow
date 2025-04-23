import psi4
 
psi4.set_memory('4000 mb')
psi4.set_num_threads(4)
 
mol_name = psi4.geometry("""
 0 1 
 C 2.92025 1.83134 1.76256 
 O 4.01415 1.03695 1.31128 
 C 3.62182 0.02385 0.46355 
 C 4.60163 -0.78687 -0.07106 
 N 4.31959 -1.80659 -0.91323 
 C 3.01822 -1.99399 -1.20843 
 C 2.68155 -3.02912 -2.06861 
 C 1.35868 -3.28347 -2.41673 
 C 0.35235 -2.48255 -1.89461 
 C 0.65942 -1.42145 -1.02286 
 C -0.43989 -0.59810 -0.49340 
 C -1.25772 -1.01485 0.58080 
 O -1.01605 -2.25669 1.13432 
 C -0.28648 -2.11489 2.35386 
 C -2.33535 -0.24069 1.03874 
 O -3.08023 -0.78280 2.05736 
 C -4.26759 -0.10156 2.43628 
 C -2.58880 1.00668 0.45679 
 C -1.76137 1.43826 -0.59689 
 C -0.70187 0.66490 -1.08252 
 C -0.10396 1.40812 -2.14460 
 C -0.79088 2.59163 -2.27806 
 N -1.78702 2.60216 -1.34072 
 C -2.73220 3.68180 -1.16079 
 C 1.99682 -1.16728 -0.66626 
 N 2.31805 -0.15649 0.17911 
 H 2.41475 2.32521 0.92511 
 H 2.21483 1.23671 2.35368 
 H 3.32350 2.61250 2.41427 
 H 5.64852 -0.63720 0.16396 
 H 3.47298 -3.65490 -2.47456 
 H 1.11571 -4.10299 -3.08777 
 H -0.68473 -2.68285 -2.15942 
 H -0.93094 -1.74728 3.15840 
 H 0.07581 -3.10651 2.64136 
 H 0.58395 -1.45959 2.23885 
 H -4.03653 0.87063 2.88356 
 H -4.77076 -0.70535 3.19768 
 H -4.95372 0.00423 1.58935 
 H -3.39911 1.64149 0.79194 
 H 0.74639 1.10471 -2.74141 
 H -0.64491 3.41800 -2.96131 
 H -2.52035 4.47867 -1.87829 
 H -3.73905 3.29120 -1.32892 
 H -2.62946 4.06502 -0.14236 
""") 
 
psi4.set_module_options('scf', {'maxiter': 75}) 
psi4.set_options({'basis': '6-31g*', 'scf_type': 'df'})
 
energy, wfn = psi4.optimize('b3lyp', engine='geometric', return_wfn=True)  
 
optimized_geometry_xyz = wfn.molecule().to_string('xyz') 
with open('12_opt.xyz', 'w') as xyz_file:
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
 
with open('12_opt_energy_and_gap.txt', 'w') as file:
    file.write(f"Optimized energy: {energy_ev:.6f} eV\n") 
    file.write(f"HOMO: {homo_energy_ev:.6f} eV\n") 
    file.write(f"LUMO: {lumo_energy_ev:.6f} eV\n") 
    file.write(f"Energy gap (HOMO-LUMO): {energy_gap_ev:.6f} eV\n") 
 
