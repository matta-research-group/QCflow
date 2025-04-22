#!/bin/bash --login
#SBATCH -o 12_sp_c.out 
#SBATCH -e 12_sp_c.err 
#SBATCH --job-name=12_sp_c 
#SBATCH -p cpu 
#SBATCH --ntasks=10
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=4000 
#SBATCH --time=5:00:00 
 
module purge 
module load cuda/10.0.130-gcc-13.2.0 
python3 12_sp_c.py 
