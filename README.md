# auto_workflow

A HTC workflow leveraging [rdkit](https://github.com/rdkit/rdkit) and [cclib](https://github.com/cclib/cclib).

- SMILES with attachment points are combined together to form dimers and trimers
- `.com` Gaussian input files, `.sh` SLURM submission files are created.
- jobs are submitted to the HPC cluster
- output files are parsed to submit further calculations or retrieve DFT descriptors

Naming Convention:

- Torsional scan neutral → tor
- Optimisation neutral/Population analysis → pop_opt_n
- Vertical anion → ver_a
- Vertical cation → ver_c
- Optimisation anion → opt_a
- Optimisation cation → opt_c
