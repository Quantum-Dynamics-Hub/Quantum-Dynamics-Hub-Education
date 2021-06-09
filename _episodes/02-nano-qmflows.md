---
title: "nano-qmflows workflows"
date: June 15, 2021, 11:00 am - 1:00 pm EDT
---


# Tutorial
<a name="toc"></a>

1. [The single_points workflow](#single_points)
2. [The absorption_spectrum workflow](#absorption_spectrum)
3. [The distribute_derivative_couplings workflow](#derivative_couplings)

## 1. The single_points workflow
<a name="single_points"></a> [Back to TOC](#toc)

A single point calculation on the relaxed geometry of a Cd33Se33 system has been performed according to the corresponding [tutorial](https://qmflows-namd.readthedocs.io/en/latest/single_points.html) in the nano-qmflows documentation.
Copy the resulting `Cd33Se33.hdf5` file in your working directory and use it to:
1. Calculate the HOMO-LUMO gap in eV.
2. Plot the energy of the Kohn-Sham orbitals considered in the active space. (Suggestion: use [matplotlib.pyplot.barh](https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.barh.html))


## 2. The absorption_spectrum workflow
<a name="absorption_spectrum"></a> [Back to TOC](#toc)

Calculate the oscillator strength of the lowest lying excited states of our Cd33Se33 system within the single orbital transitions approximation.
To do that, in your working directory:
- copy the pre-computed `Cd33Se33.hdf5` file;
- copy the file containing the coordinates of relaxed Cd33Se33 geometry in an xyz format, `Cd33Se33.xyz`;
- create an `absorption_spectrum_Cd33Se33.yml` input file and customize it according to the previous requirements using to this [tutorial](https://qmflows-namd.readthedocs.io/en/latest/absorption_spectrum.html);
- create a `launch.sh` submission script containing:

      #!/bin/sh
      #SBATCH --partition=valhalla  --qos=valhalla
      #SBATCH --clusters=faculty
      #SBATCH --time=00:10:00
      #SBATCH --nodes=1
      #SBATCH --ntasks-per-node=12
      #SBATCH --mem=12000
       
      eval "$(/projects/academic/cyberwksp21/Software/Conda/Miniconda3/bin/conda shell.bash hook)"
      conda activate qmflows
      module load cp2k/8.1-sse
       
      run_workflow.py -i absorption_spectrum_Cd33Se33.yml
       
and submit your calculation with:
 
    sbatch launch.sh
    
1. How many singly excited configurations do you expect to be calculated?
2. Use the script [convolution.py](https://github.com/SCM-NV/nano-qmflows/blob/master/scripts/qmflows/convolution.py) to plot the corresponding absorption spectrum.

## 3. The distribute_derivative_couplings workflow
<a name="#derivative_couplings"></a> [Back to TOC](#toc)

Follow the [Derivative_Couplings tutorial](https://qmflows-namd.readthedocs.io/en/latest/derivative_couplings.html) to calculate the derivative couplings of 
