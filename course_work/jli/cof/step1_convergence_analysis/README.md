CP2K convergence analysis for subNEUCOF1

The convegence analsysis benchmarks the cut off energy (CUTOFF),
the number of integration grid points (NGRIDS), and the relative 
cut off energy (REL_CUTOFF) for the single point caculation.

The method uses the PBE functional and the basis sets and potential are:

Kind    Basis set        Potential
 H   DZVP-MOLOPT-GTH    GTH-PBE-q1
 C   DZVP-MOLOPT-GTH    GTH-PBE-q4
 N   DZVP-MOLOPT-GTH    GTH-PBE-q5
 B   DZVP-MOLOPT-SR-GTH GTH-PBE-q3
 S   DZVP-MOLOPT-GTH    GTH-PBE-q6

The scan values of CUTOFF, NGRIDS, and REL_CUTOFF are:
CUTOFF:      50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1200
NGRIDS:      4,5,6,7,8,9,10,11,12,13,14
REL_CUTOFF:  50,60,70,80,90,100,110,120,130,140,150

Each subfolder contains:
1. cof.inp and cof.xyz, which are the input file and coordinates used for the calculation
2. a basis set and potential file used for the calculation
3. run_cp2k_convergence_analysis.sh, which is a shell script for convergence analysis
4. check_results.sh, which is a shell script for results collection
5. energy.txt, which saves the 1-step SCF energy
6. summary.txt, which saves the detials, e.g., scan variables, grids point partitions, timing.
7. log files, named as cof-${ngrids}-${cutoff}-${rel_cutoff}.log
