#!/bin/bash -l
#SBATCH --job-name=asymp

# Specify the name and location of i/o files.  
# “%j” places the job number in the name of those files.
#SBATCH --output=z_err_out/asymp_%j.out
#SBATCH --error=z_err_out/asymp_%j.err

# Send email notifications.  
# #SBATCH --mail-type=END # other options are ALL, NONE, BEGIN, FAIL
# #SBATCH --mail-user=aenoble@ucdavis.edu

# Specify the partition.
#SBATCH --partition=hi # other options are low, med, bigmem, serial.

# Specify the number of requested nodes.
#SBATCH --nodes=1
#SBATCH --exclude=c9-70

# Specify the number of tasks.
#SBATCH --ntasks=1 

# Specify the number of jobs to be run, 
# each indexed by an integer taken from the interval given by “array”.
#SBATCH --array=0-215

hostname
echo "SLURM_NODELIST = $SLURM_NODELIST"
echo "SLURM_NODE_ALIASES = $SLURM_NODE_ALIASES"
echo "SLURM_NNODES = $SLURM_NNODES"
echo "SLURM_TASKS_PER_NODE = $SLURM_TASKS_PER_NODE"
echo "SLURM_NTASKS = $SLURM_NTASKS"
echo "SLURM_JOB_ID = $SLURM_JOB_ID"
echo "SLURM_ARRAY_TASK_ID = $SLURM_ARRAY_TASK_ID"

./a.out $SLURM_ARRAY_TASK_ID

