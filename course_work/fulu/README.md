The files for my projects using the packages covered in this workshop are listed in this directory.

I mainly adopted two packages, i.e., NEXMD and Libra, to conduct my non-adiabatic excited state dynamics (NAESMD) simulations. For the study using Libra, the Quantum Espresso has been used to perform MD simulations and single point calculations.

As included in "1.nexmd_project", I have used NEXMD to simulate the NAESMD in two building units of 1D molecular wires.

As included in "2.libra_QE_project", I have used Libra-QE to study the NAESMD in a monomer and a dimer of ZnTHPP ((5,10,15,20-Tetrakis(4-hydroxyphenyl)-zinc-porphine) molecule. The jupyter notebooks used in this project was copied from the tutorial of Libra (https://github.com/compchem-cybertraining/Tutorials_Libra/tree/master/6_dynamics/2_nbra_workflows).

Please note that the trajectories from NEXMD and the single-point calculation results from QE are not included in the repository due to their huge sizes. Those files are available on request by contact me via fzheng@uni-bremen.de.

The NEXMD calculations are performed on the clusters of Max Planck Institute for the Physics of Complex Systems (MPI-PKS) and the calculations with Libra and QE are performed on the Draco cluster of Max Planck Computing and Data Facility (MPCDF). Computational resources and technique support from MPI-PKS and MPCDF are gratefully acknowledged.
