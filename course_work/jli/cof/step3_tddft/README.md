CP2K TDDFT calculation for subNEUCOF1

The method uses the PBE functional and the basis set and potentials are:
Kind    Basis set        Potential
 H   DZVP-MOLOPT-GTH    GTH-PBE-q1
 C   DZVP-MOLOPT-GTH    GTH-PBE-q4
 N   DZVP-MOLOPT-GTH    GTH-PBE-q5
 B   DZVP-MOLOPT-SR-GTH GTH-PBE-q3
 S   DZVP-MOLOPT-GTH    GTH-PBE-q6

This folder contains
1. an input file used for the calculation
2. a basis set and potential file used for the calculation
3. a shell script for slurm job submission
4. COF-MOS-1_0.molden, molecular orbital file

Note the DFT calculation is done but the TD-DFT calculation is
terminated due to out-of-memory. But this is a good exercise to
run the TD-DFT calculation.
