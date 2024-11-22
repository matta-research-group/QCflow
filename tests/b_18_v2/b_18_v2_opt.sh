#!/bin/bash --login
#SBATCH -o b_18_v2_opt.out 
#SBATCH -e b_18_v2_opt.err 
#SBATCH --job-name=b_18_v2_opt 
#SBATCH -p cpu 
#SBATCH --ntasks=20
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=4000 
#SBATCH --time=48:00:00 
 
INPUTFILE=b_18_v2_opt.com 
OUTPUTFILE=b_18_v2_opt.log 
 
module purge 
module load gaussian_sse4/16-C-gcc-13.2.0 
export GOMP_CPU_AFFINITY=$SGE_BINDING 
export KMP_AFFINITY="explicit,proclist=$SGE_BINDING,verbose" 
#source $g16root/bsd/g16.login 
 
echo "G16 job \$SLURM_JOBID" 
echo "INPUT \$INPUTFILE" 
echo "OUTPUT \$OUTPUTFILE" 
echo "Running \$SLURM_NTASKS on \$SLURM_JOB_NODELIST" 
 
#Execution Line 
g16 $INPUTFILE > $OUTPUTFILE 
