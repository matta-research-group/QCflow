#!/bin/bash
#SBATCH --job-name=QCflow-env
#SBATCH -o QCflow-env.out 
#SBATCH -e QCflow-env.err 
#SBATCH -p cpu
#SBATCH --ntasks=32
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=8000 
#SBATCH --time=24:00:00
#SBATCH --mail-user=k2255489@kcl.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL    

MICROMAMBA=/cephfs/volumes/hpc_data_usr/k2255489/286fa9ad-d272-4b7d-a3ee-ccc6df026c4c/bin/micromamba
export MAMBA_ROOT_PREFIX=/cephfs/volumes/hpc_home/k2255489/6adaf325-0832-46f3-9840-3cdcb52285e8

eval "$($MICROMAMBA shell hook --shell=bash)"

$MICROMAMBA env create -f QCflow.yml -y
