#!/bin/bash --login
#SBATCH -o 12_anion.out 
#SBATCH -e 12_anion.err 
#SBATCH --job-name=12_anion 
#SBATCH -p cpu 
#SBATCH --ntasks=10
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=4000 
#SBATCH --time=5:00:00 
 
module purge 
module load cuda/10.0.130-gcc-13.2.0 
python3 12_anion.py 
