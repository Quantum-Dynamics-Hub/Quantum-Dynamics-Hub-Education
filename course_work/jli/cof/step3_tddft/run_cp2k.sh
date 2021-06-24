#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1
#SBATCH --time=21:00:00
#SBATCH --job-name=cp2k-t
#SBATCH --partition=valhalla  --qos=valhalla
#SBATCH --clusters=faculty
#SBATCH --mem=200GB
#SBATCH --output=%j.o.slurm
#SBATCH --error=%j.e.slurm


export INPUT=cof
export WORKDIR=/panasas/scratch/grp-cyberwksp21/ub2053/cof/tddft

module load cp2k/8.1-sse
#module load openmpi/3.0.3/gcc-7.3.0

export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

export OMP_NUM_THREADS=1
export NPROCS=$SLURM_NTASKS

cd $WORKDIR
> $INPUT.log
mpirun -n $NPROCS cp2k.popt -i $INPUT.inp -o $INPUT.log

    

