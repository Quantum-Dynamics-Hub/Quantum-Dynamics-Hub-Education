#!/bin/bash -l
#SBATCH -J ZnTHPP_relax            # job name
#SBATCH -o ./job.out.%j      # standard out file
#SBATCH -e ./job.err.%j      # standard err file
#SBATCH -D ./                # work directory
#SBATCH --partition=general
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu 1024
#SBATCH --time=24:00:00       # time limit, max 24h
module load qe/6.4
srun pw.x -inp ZnTHPP_dimer_scf.in > ZnTHPP_dimer_scf.out
