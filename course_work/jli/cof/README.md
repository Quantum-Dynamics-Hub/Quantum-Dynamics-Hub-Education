Author: Jingbai Li@Northeastern University

Project description

In 2020, the Lopez Lab designed a new subphthalocyanine-based 3D 
covalent organic framework (NEUCOF1) capable of forming cocrystals 
with fullerene (C60) via periodic ball-and-socket binding motifs 
(J. Phys. Chem. C 2020, 124, 9126−9133). NEUCOF1 has a high 
cocrystalline surface area and long-range order that eliminates 
the typical surface area vs. long-range order trade-off in organic 
photovoltaics (OPVs). The plane-wave density functional theory 
calculation suggests exciton charge transfer from NEUCOF1 to the 
pocket-bound fullerenes, followed by a subsequent free electron 
transfer to the nanowire of C60 acceptors. 

In this workshop, I plan to use CP2K to explore the hot electron 
relaxation process through the NECUCOF1-C60 donor-acceptor interface.
The model creates a subsystem of NEUCOF1 (subNEUCOF1) that contains 
two COF shell enclosing a C60 molecule. 

The project includes the followig steps. This workshop project only covers
step 1-3 due to the overly buzy excercise time.

1. convergence analysis
2. geometry optimization
3. time-dependent density functional theory calculation
4. excited-state dynamics of NEUCOF1
5. interfacing CP2K to machine-learning nonadiabatic molecular dyanmics
   code - PyRAI2MD(https://github.com/lopez-lab/PyRAI2MD)
