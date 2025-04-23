#!/bin/bash --login
#SBATCH -o 1_opt.out 
#SBATCH -e 1_opt.err 
#SBATCH --job-name=1_opt 
#SBATCH -p cpu 
#SBATCH --ntasks=10
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=4000 
#SBATCH --time=24:00:00 
 
module purge 
module load cuda/10.0.130-gcc-13.2.0 
python3 1_opt.py 
