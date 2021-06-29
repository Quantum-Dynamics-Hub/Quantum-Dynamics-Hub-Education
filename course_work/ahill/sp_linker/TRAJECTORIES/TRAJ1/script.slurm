#!/bin/sh
#SBATCH --partition=valhalla --qos=valhalla
#SBATCH --clusters=faculty
#SBATCH --account=cyberwksp21
#SBATCH --time=10:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=1500

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST="$SLURM_JOB_NODELIST
echo "SLURM_NNODES="$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "working directory="$SLURM_SUBMIT_DIR

eval "$(/projects/academic/cyberwksp21/Software/Conda/Miniconda3/bin/conda shell.bash hook)"
conda activate libra

NPROCS=`srun --nodes=${SLURM_NNODES} bash -c 'hostname' |wc -l`
echo NPROCS=$NPROCS
export I_MPI_PMI_LIBRARY=/usr/lib64/libmpi.so

export OMP_NUM_THREADS=8
export NX=/projects/academic/cyberwksp21/Software/NX-2.2-B01/bin
export DFTBP=/projects/academic/cyberwksp21/Software/dftbplus-17.1.x86_64-linux/bin
export TDDFTB=/projects/academic/cyberwksp21/Software/dftbplus-17.1.x86_64-linux/bin
export DP_TOOLS=/projects/academic/cyberwksp21/Software/dftbplus-17.1/tools/dptools/bin
export PATH=$DP_TOOLS/:$PATH
export PYTHONPATH=$PYTHONPATH:/projects/academic/cyberwksp21/Software/dftbplus-17.1/tools/dptools/lib/python2.7/site-packages/
export PATH=$PATH:/projects/academic/cyberwksp21/Software/dftbplus-17.1/tools/dptools/bin
$NX/moldyn.pl > moldyn.log 
