# auto_workflow

A HTC workflow leveraging [rdkit](https://github.com/rdkit/rdkit) and [cclib](https://github.com/cclib/cclib).

- SMILES with attachment points are combined together to form dimers and trimers
- `.com` Gaussian input files, `.sh` SLURM submission files are created.
- jobs are submitted to the HPC cluster
- output files are parsed to submit further calculations or retrieve DFT descriptors
