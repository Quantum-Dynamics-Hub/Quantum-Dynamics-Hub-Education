#!/bin/sh
#SBATCH --account=cyberwksp21
#SBATCH --partition=valhalla  --qos=valhalla
#SBATCH --clusters=faculty

#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=12000

module load columbus

runc -m 12000 > runls
