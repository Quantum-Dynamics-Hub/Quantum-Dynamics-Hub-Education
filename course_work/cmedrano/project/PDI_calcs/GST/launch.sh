#!/bin/sh
#SBATCH --partition=valhalla  --qos=valhalla
#SBATCH --account=cyberwksp21

#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=12000

module load nexmd

nexmd.exe > output
