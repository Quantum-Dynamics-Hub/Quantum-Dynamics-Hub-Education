#!/bin/sh
#SBATCH --partition=valhalla  --qos=valhalla
#SBATCH --clusters=faculty
#SBATCH --account=cyberwksp21

#SBATCH --time=00:45:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=12000

eval "$(/projects/academic/cyberwksp21/Software/Conda/Miniconda3/bin/conda shell.bash hook)"
conda activate qmflows
module load cp2k/8.1-sse

init_armc fitting.yaml
