#!/bin/sh
#SBATCH --partition=valhalla  --qos=valhalla
#SBATCH --clusters=faculty
#SBATCH --account=cyberwksp21

#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=12000
 
eval "$(/projects/academic/cyberwksp21/Software/Conda/Miniconda3/bin/conda shell.bash hook)"
conda activate qmflows
module load cp2k/8.1-sse
 
run_workflow.py -i absorption_spectrum_Cd33Se33.yml
