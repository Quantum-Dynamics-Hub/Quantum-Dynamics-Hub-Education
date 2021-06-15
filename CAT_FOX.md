---
title: "CAT workflow"
date: June 17, 2021, 11:00 am - 1:00 pm EDT
---


# Tutorial
<a name="toc"></a>

0. [Setup](#setup)
1. [The qd_build workflow](#qd_build)
2. [The fitting workflow](#fitting)

## 0. Setup
<a name="setup"></a> [Back to TOC](#toc)

In your working directory, copy the folder containing all the files required for this tutorial:

    cp -r /projects/academic/cyberwksp21/Instructors_material/rpascazio/exercises/ .

Please refer to the following documentations:
1. The Building a Quantum Dot Model assignment is based on the [CAT documentation](https://cat.readthedocs.io/en/latest/) and on the [Building a Quantum Dot Model Tutorial](https://nanotutorials.readthedocs.io/en/latest/1_build_qd.html);
2. The Forcefield Optimization assignment is based on the [Auto-FOX documentation](https://auto-fox.readthedocs.io/en/latest/) and on the [Forcefield Optimization Tutorial](https://nanotutorials.readthedocs.io/en/latest/2_fitting.html).

## 1. Building a Quantum Dot Model
<a name="qd_build"></a> [Back to TOC](#toc)

We aim to create nanocrystals from a charge-balanced Cd68Se68 core. Use the corresponding `CdSe.xyz` file in your `1_qd_build` directory to:
1. Replace a 30% fraction of the Se ions in the model with Cl dummies in a file called `CdSe_30Cl.xyz` with random distribution and use the file to replace the Cl on the surface with stearate anions.
2. Replace a 20% fraction of the Cd ions in the model with Na dummies in a file called `CdSe_20Na.xyz` with clustered distribution and use the file to replace the Na on the surface with oleylammonium cations.
(Suggestion: their SMILES strings are available on [PUBCHEM](https://pubchem.ncbi.nlm.nih.gov/)).

## 2. The fitting workflow
<a name="fitting"></a> [Back to TOC](#toc)

We aim to fit the classical forcefield parameters of a CsPbBr3 core from a previous QM-MD trajectory using the Adaptive Rate Monte Carlo (ARMC) algorithm.
1. Modify the yaml script to fit the `CsPbBr3_MD.xyz` trajectory (you can find them in your `2_fitting` directory).
Use the oxidation states of the atoms (respectively 1.0, 2.0, -1.0 for Cs, Pb, Br) as starting parameters for their charges, and modify the `guess_rdf.py` script in your `scripts` directory to obtain the starting values for the sigmas. (Beware of the units!). Run the fitting.
2. Run the parametrization for around 5 steps (Suggestion: use the  `logfile` to check the number of last iteration). Use the python scripts in the `scripts` subdirectory to plot the errors of the ARMC and the Radial Distribution Functions (RDFs) of the best step from the `armc.hdf5` file. Repeat the same procedure after 15 more steps (around 20 steps in total). How do you expect the plots to change over time?
