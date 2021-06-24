#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1
#SBATCH --time=21:00:00
#SBATCH --job-name=cp2k-t
#SBATCH --partition=valhalla  --qos=valhalla
#SBATCH --clusters=faculty
#SBATCH --mem=100GB
#SBATCH --output=%j.o.slurm
#SBATCH --error=%j.e.slurm


export INPUT=cof
export WORKDIR=/panasas/scratch/grp-cyberwksp21/ub2053/cof/convergence_analysis/XYZ/relcutoff

module load cp2k/8.1-sse
#module load openmpi/3.0.3/gcc-7.3.0

export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

export OMP_NUM_THREADS=1
export NPROCS=$SLURM_NTASKS

cd $WORKDIR

> energy.txt
#for ngrids in {4,6,8,10,12,14}
#do
    ngrids=6
    #for cutoff in {50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800} 
#    for cutoff in {850,900,950,1000,1050,1100,1150,1200}
#    do
    cutoff=1000
        for relcutoff in {50,60,70,80,90,100,110,120,130,140,150}
        do
           # Substituting the CUTOFF value 
           sed -i "/\<CUTOFF\> /c\     CUTOFF $cutoff" $INPUT.inp
           # Substituting the CUTOFF value
           sed -i "/\<REL_CUTOFF\> /c\    REL_CUTOFF $relcutoff" $INPUT.inp
           # Substituting the NGRIDS value 
           sed -i "/\<NGRIDS\> /c\     NGRIDS $ngrids" $INPUT.inp

           > $INPUT.log
           mpirun -n $NPROCS cp2k.popt -i $INPUT.inp -o $INPUT.log
           scf=`grep -A4 'SCF WAVEFUNC' $INPUT.log|tail -n1|awk '{print $6}'`

           echo $ngrids $cutoff $relcutoff $scf >> energy.txt
           mv $INPUT.log $INPUT-$ngrids-$cutoff-$relcutoff.log
        done
#    done
#done



    

