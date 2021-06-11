# Scripts for Libra/qmflows setup

These scripts were used to set up the Libra installation, conda
environments, and the libra-plus and qmflows Jupyter Notebooks kernels at CCR.

`libra_env_recipe.sh`
    
   Installation procedure for libra prerequisites and conda environment
   Refer to libra documentation for further information.

`qmflows_env_recipe.sh`
   
   Installation procedure for the qmflows conda environment.

`prepend_and_launch.sh`
    
   These scripts are called upon launch of the Jupyter libra-plus
   and qmflows kernels. They set LD_LIBRARY_PATH, PATH, and various 
   environment variables to enable Libra to be used with various compiled
   codes on the compute cluster.
   
   The prepend_and_launch.sh scripts should be placed in the kernel
   directory for the respective kernel, and specified in the arguments of
   the appropriate kernel.json so it will be called.
   
`kernel.json`

    Example kernel.json for Libra, showing specification of launch script.       

## Compiled codes

The following codes have been installed on CCR computing resources in 
support of this project:

 - columbus
 - cp2k
 - dftb+
 - dynemol
 - eQE
 - ergoSCF (serial and MPI)
 - LAMMPS
 - MOPAC
 - Newton-X
 - Q-Chem
 - Quantum Espresso 6.2.1
 - QXMD

## Compilers
 
The following compilers were used for all codes listed above, as appropriate: 
 
 - MKL 2020.2
 - INTEL 18.3
 - openmpi 3.0.3/gcc 7.3.0

## Computing resources

Codes listed above were compiled specifically for CCR's faculty cluster `valhalla` partition.
Node specs (two kinds) are as follows:

CPU-E5-2620v3:

- CPUs 12
- INTEL
- Memory 128000 (MB) 

CPU-E5-2650v4:

- CPUs 24
- INTEL
- Memory 256000 (MB)
 


