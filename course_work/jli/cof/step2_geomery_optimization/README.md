CP2K geometry optimization for subNEUCOF1

The optimization only minimize the nuclear position in the unit cell.
The cell lattice vectors are obtained from the previous optimized NEUCOF1
(J. Phys. Chem. C 2020, 124, 9126−9133)

The lattice vectors are:
    A  23.01504              0                0              
    B -12.3164073652008     20.1158047946615  0              
    C  -0.363649935235012    0.26423684160044 38.8413789644053

Note the subNEUCOF1 is a 2D system. The C vector was doubled to ensure 
adequate vacuum space along the z-axis.

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
3. dftd3.dat for Grimme's D3 correction
4. a shell script for slurm job submission
5. a log file and optimized geometry
6. COF-frc-1.xyz, the latest geometry
