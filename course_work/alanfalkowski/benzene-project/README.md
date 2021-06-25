----------------------------------------------
Writed by Alan Guilherme Falkowski.

email: agf18@ifi.unicamp.br
----------------------------------------------

Description:
-----------

Here we have the inputs, outputs and so on.

In general, all calculations has been performed with the aug-cc-pVDZ basis set and the point group used is the D2h.


Acronyms used:
-------------

HF: Hartree-Fock

SCF: Self Consistent Field

MP2: Møller-Plesset perturbation theory of second order

CAS: Complete Active Space

MRCI: Multi Reference Configutation Interaction

MRCISD: MRCI with Singles and Doubles excitations

EOM-CCSD: Equation-of-motion Coupled Cluster with Singles and Doubles excitations


CONTENTS:
----------------------------------------------

1) C6H6_GAMESS_MP2_OPT_d2h: Geometry optimization performed in the GAMESS-US package, with the MP2 approach.

2) C6H6_GAMESS_HF_SCF_d2h: Hartree-Fock SCF calculation. This calculations helped me to know about the occupations of each irreducible representation of D2h point group.

3) C6H6_CAS_6-6_d2h_SP: Single point calculation with the CAS-SCF approach, performed in the COLUMBUS package. The active space used is CAS(6,6), i.e., 6 electrons in 6 orbitals. I used all pi occupied orbitals (3) and the first 3 pi unnocupied orbitals.

4) C6H6_MRCISD_6-6_d2h_SP: Single point calculation with the MR-CISD method, performed in the COLUMBUS package. The active space used is MRCI(6,6), the same electrons and orbitals used in the CAS(6,6). More details about the table of orbitals are in the file "C6H6_CAS_6-6_d2h_table.odt". 

5) C6H6_PSI4_EOMCCSD_d2h: Also I performed an EOM-CCSD calculation to compare with my multireference calculations. This has been performed using the psi4 package. 

6) C6H6_CAS_6-6_d2h_table.pdf: Informations about occupations, active space, etc.

7) detailed_options_columbus.pdf: A step-by-step of how I have done the calculations in COLUMBUS, for both CASSCF and MRCI.

8) launch.sh: An example SLURM submission file for clusters that use SLURM

