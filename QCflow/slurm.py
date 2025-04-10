import subprocess

def write_slurm(job_name, mol_name, cpus=10):
    """
    Writes a SLURM batch script for a specified job type and molecule name.

    Parameters
    ----------
    job_name (str): The type of job to run. Possible values:
        - 'sp': Single Point neutral
        - 'opt': Optimisation neutral
        - 'tor': Torsional scan neutral                                                                                
        - 'pop_opt_n': Optimisation neutral + Population analysis                                                                                
        - 'sp_a': Single point anion                                                                                
        - 'sp_c': Single point cation
        - 'opt_a': Optimisation anion
        - 'opt_c': Optimisation cation
        - 'n_a_geo': Neutral charge, optimised anion geometry
        - 'n_c_geo': Neutral charge, optimised cation geometry
        - 'sp_hirsh': Single Point Hirshfeld 
    mol_name (str): The name of the dimer from the dictionary, e.g., if fragment 0 was attached to fragment 1,
                    then the dimer name is '0_1'.
    cpus (int, optional): The number of CPUs to allocate for the job. Default is 10.
    
    Notes
    -----
    The function generates a SLURM batch script file named '{mol_name}_{job_name}.sh' with appropriate
    configurations based on the job type and molecule name. The script includes settings for job name,
    output and error files, partition, number of tasks, nodes, CPUs per task, memory per CPU, and time limit.
    It also sets up the environment and execution line for running Gaussian 16 (g16) with the specified input
    and output files.
    """

    file_name = f'{mol_name}_{job_name}.sh'
    
    if (job_name == 'sp') or (job_name == 'ver_a') or (job_name == 'ver_c') or (job_name == 'sp_hirsh') or (job_name == 'n_a_geo') or (job_name == 'n_c_geo'):
        calc_time = '24:00:00'
    else:
        calc_time = '48:00:00'

    title = f'#!/bin/bash --login'
    with open(file_name, 'w') as file:
        file.write(f'{title}\n')#
        file.write(f'#SBATCH -o {mol_name}_{job_name}.out \n')
        file.write(f'#SBATCH -e {mol_name}_{job_name}.err \n')#
        file.write(f'#SBATCH --job-name={mol_name}_{job_name} \n')
        file.write(f'#SBATCH -p cpu \n')
        file.write(f'#SBATCH --ntasks={cpus}\n')
        file.write(f'#SBATCH --nodes=1 \n')
        file.write(f'#SBATCH --cpus-per-task=1 \n')
        file.write(f'#SBATCH --mem-per-cpu=4000 \n')
        file.write(f'#SBATCH --time={calc_time} \n') #reduced to speed up queue time
        file.write(' \n')#
        file.write(f'INPUTFILE={mol_name}_{job_name}.com \n')
        file.write(f'OUTPUTFILE={mol_name}_{job_name}.log \n')
        file.write(' \n')
        file.write(f'module purge \n')
        file.write(f'module load gaussian_sse4/16-C-gcc-13.2.0 \n')
        file.write(f'export GOMP_CPU_AFFINITY=$SGE_BINDING \n')
        file.write(f'export KMP_AFFINITY="explicit,proclist=$SGE_BINDING,verbose" \n')
        file.write(f'#source $g16root/bsd/g16.login \n')
        file.write(' \n')
        file.write(f'echo "G16 job \$SLURM_JOBID" \n')
        file.write(f'echo "INPUT \$INPUTFILE" \n')
        file.write(f'echo "OUTPUT \$OUTPUTFILE" \n')
        file.write(f'echo "Running \$SLURM_NTASKS on \$SLURM_JOB_NODELIST" \n')
        file.write(' \n')
        file.write(f'#Execution Line \n')
        file.write(f'g16 $INPUTFILE > $OUTPUTFILE \n')

def write_slurm_psi4(job_name, mol_name, time=24, cpus=10):
    """
    Writes a SLURM batch script for a specified job type and molecule name.

    Parameters
    ----------
    job_name (str): The type of job to run. Possible values:
        - 'sp': Single Point neutral
        - 'opt': Optimisation neutral
        - 'cation': Geometry optimisation cation (opt_c) and single of neutral charge, cation geometry (n_c_geo)
        - 'anion': Geometry optimisation anion (opt_a) and single of neutral charge, anion geometry (n_a_geo)
        - 'sp_c': Single point calculation of neutral geometry at cation charge
        - 'sp_a': Single point calculation of neutral geometry at anion charge
    mol_name (str): The name of the dimer from the dictionary, e.g., if fragment 0 was attached to fragment 1,
                    then the dimer name is '0_1'.
    time (int, optional): The time limit for the job in hours. Default is 24. (Max is 48)
    cpus (int, optional): The number of CPUs to allocate for the job. Default is 10.
    
    Notes
    -----
    The function generates a SLURM batch script file named '{mol_name}_{job_name}.sh' with appropriate
    configurations based on the job type and molecule name. The script includes settings for job name,
    output and error files, partition, number of tasks, nodes, CPUs per task, memory per CPU, and time limit.
    It also sets up the environment and execution line for running Psi4 job with the specified input
    and output files.
    """

    file_name = f'{mol_name}_{job_name}.sh'
    
    calc_time = f'{time}:00:00'

    title = f'#!/bin/bash --login'
    with open(file_name, 'w') as file:
        file.write(f'{title}\n')#
        file.write(f'#SBATCH -o {mol_name}_{job_name}.out \n')
        file.write(f'#SBATCH -e {mol_name}_{job_name}.err \n')#
        file.write(f'#SBATCH --job-name={mol_name}_{job_name} \n')
        file.write(f'#SBATCH -p cpu \n')
        file.write(f'#SBATCH --ntasks={cpus}\n')
        file.write(f'#SBATCH --nodes=1 \n')
        file.write(f'#SBATCH --cpus-per-task=1 \n')
        file.write(f'#SBATCH --mem-per-cpu=4000 \n')
        file.write(f'#SBATCH --time={calc_time} \n') #reduced to speed up queue time
        file.write(' \n')
        file.write(f'module purge \n')
        file.write(f'module load cuda/10.0.130-gcc-13.2.0 \n')
        file.write(f'python3 {mol_name}_{job_name}.py \n')


def is_job_in_queue(submission_name):
    """
    Checks if the job has been successfully submitted into the queue

    Parameters
    ----------
    submission_name (str): The name of the job that has been submitted to the HPC

    Returns
    -------
    True or False
    """
    try:
        result = subprocess.run(
            ['squeue', '--me'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        # Search for job_name in the output
        for line in result.stdout.splitlines():
            if submission_name in line:
                return True
        return False
    except subprocess.CalledProcessError as e:
        print(f"Error checking squeue: {e.stderr.strip()}")
        return False


def submit_slurm_job(job_name, mol_name, max_retries=5, wait_seconds=30):
    """
    Submits a SLURM job using the specified job name and molecule name. Works on the KCL CREATE HPC.

    Parameters
    ----------
    job_name (str): The type of job to run. Possible values:
        - 'sp': Single Point neutral
        - 'opt': Optimisation neutral
        - 'tor': Torsional scan neutral                                                                                
        - 'pop_opt_n': Optimisation neutral + Population analysis                                                                                
        - 'sp_a': Single point anion                                                                                
        - 'sp_c': Single point cation
        - 'opt_a': Optimisation anion
        - 'opt_c': Optimisation cation
        - 'n_a_geo': Neutral charge, optimised anion geometry
        - 'n_c_geo': Neutral charge, optimised cation geometry
        - 'sp_hirsh': Single Point Hirshfeld 

    mol_name (str): The name of the dimer from the dictionary. For example, if fragment 0 was attached to fragment 1,
                    then the dimer name would be '0_1'.

    max_retries (int): The maximum amount of times a job will attempt to submit. Deafult is 5.

    wait_seconds (int): How long python will go to sleep inbetween attempts to submit

    Returns
    -------
    bytes: The standard output from the SLURM job submission command.
    """

    string = f'sbatch {mol_name}_{job_name}.sh'

    submission_name = f'{mol_name}_{job_name}.sh'

    for attempt in range(1, max_retries + 1):
        if is_job_in_queue(submission_name):
            print(f"Job {submission_name} is already in the queue. Skipping submission.")
            return None

        try:
            process = subprocess.run(
                string,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                shell=True,
                check=True
            )
            print(f"Successfully submitted job {mol_name}_{job_name}")
            return process.stdout

        except subprocess.CalledProcessError as e:
            print(f"[Attempt {attempt}] Error submitting job {mol_name}_{job_name}: {e.stderr.decode().strip()}")
            if attempt < max_retries:
                print(f"Waiting {wait_seconds} seconds before checking and retrying...")
                time.sleep(wait_seconds)
            else:
                print("Max retries reached. Moving on.")
                return None
