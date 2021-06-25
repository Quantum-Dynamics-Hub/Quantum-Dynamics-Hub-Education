#!/bin/bash

################
#
# Setting slurm options
#
################

# lines starting with "#SBATCH" define your jobs parameters

# requesting the type of node on which to run job
#SBATCH --partition long #medium #

# telling slurm how many instances of this job to spawn (typically 1)
#SBATCH --ntasks 1

# setting number of CPUs per task (1 for serial jobs)
#SBATCH --cpus-per-task 1

# export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"

# setting memory requirements
#SBATCH --mem-per-cpu 1024

# propagating max time for job to run
##SBATCH --time <days-hours:minute:seconds>
##SBATCH --time = 47:59:59
#SBATCH --time 20160 #2880 #120

# Setting the name for the job
#SBATCH --job-name ZnPc

# setting notifications for job
# accepted values are ALL, BEGIN, END, FAIL, REQUEUE
#SBATCH --mail-type  FAIL

# telling slurm where to write output and error
#SBATCH --output %j.out
#SBATCH --error %j.out

################
#
# copying your data to /scratch
#
################

# create local folder on ComputeNode
scratch=/scratch/$USER/$SLURM_JOB_ID
mkdir -p $scratch

# copy all your NEEDED data to ComputeNode
dir="$(pwd)"
echo "We are in $(pwd)"
cp $(pwd)/* $scratch
cd $scratch

# dont access /home after this line

# if needed load modules here
module add python3

# if needed add export variables here

################
#
# run the program
#
################
/data/finite/fuluzheng/6.NEXMD/NEXMD-cleaned/nexmd.exe > md.out

# copy results to data
cp -r . $dir/
echo "files copied back to cwd!"
cd

# clean up scratch
rm -rf $scratch
unset scratch

exit 0
