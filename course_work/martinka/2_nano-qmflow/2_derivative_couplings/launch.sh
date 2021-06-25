#!/bin/sh
#SBATCH -p valhalla --qos=valhalla
#SBATCH --account=cyberwksp21
#SBATCH --clusters=faculty
#SBATCH -t 00:10:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=12000
eval "$(/projects/academic/cyberwksp21/Software/Conda/Miniconda3/bin/conda shell.bash hook)"
conda activate qmflows
module load cp2k/8.1-sse

run_workflow.py -i input.yml

mkdir -p results_chunk_0
cp -r /panasas/scratch/grp-cyberwksp21/ub2047/demo/distribute_derivative_couplings_ethylene/scratch_chunk_0/chunk_0.hdf5 /panasas/scratch/grp-cyberwksp21/ub2047/demo/distribute_derivative_couplings_ethylene/scratch_chunk_0/hamiltonians/Ham_{0..-11}_* results_chunk_0
