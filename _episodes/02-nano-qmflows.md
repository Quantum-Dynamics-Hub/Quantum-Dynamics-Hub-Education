---
title: "nano-qmflows workflows"
date: June 15, 2021, 11:00 am - 1:00 pm EDT
---


# Tutorial
<a name="toc"></a>
0. [Setup](#setup)
1. [The single_points workflow](#single_points)
2. [The absorption_spectrum workflow](#absorption_spectrum)
3. [The distribute_derivative_couplings workflow](#derivative_couplings)

## 0. Setup
<a name="setup"></a> [Back to TOC](#toc)

In your working directory, copy the folder containing all the files required for the following assignments:

`cp -r /projects/academic/cyberwksp21/Instructors_material/jzito/nano-qmflows/`

Please refer to the [nano-qmflows’s documentation](https://qmflows-namd.readthedocs.io/en/latest/).


## 1. The single_points workflow
<a name="single_points"></a> [Back to TOC](#toc)

A single point calculation on the relaxed geometry of a Cd33Se33 system has been performed according to the [Single points calculation's tutorial](https://qmflows-namd.readthedocs.io/en/latest/single_points.html) (see the corresponding input file in your `1_single_points` directory).
Use the provided `Cd33Se33.hdf5` file to:
1. Calculate the HOMO-LUMO gap in eV.
2. Plot the energy (in eV) of the Kohn-Sham orbitals considered in the active space. (Suggestion: use [matplotlib.pyplot.barh](https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.barh.html))


## 2. The absorption_spectrum workflow
<a name="absorption_spectrum"></a> [Back to TOC](#toc)

Calculate the oscillator strength of the lowest lying excited states of our Cd33Se33 system within the single orbital transitions approximation.

To do that, edit the input file `absorption_spectrum_Cd33Se33.yml` provided in the directory `2_absorption_spectrum` according to the previous requirements (consult the tutorial [Absorption Spectrum](https://qmflows-namd.readthedocs.io/en/latest/absorption_spectrum.html)), then submit your calculation using the `launch.sh` submission script. 

Once the calculation is completed, copy locally the result file `output_0_sing_orb.txt` from your scratch directory and interpret it using the last part of the [tutorial](https://qmflows-namd.readthedocs.io/en/latest/absorption_spectrum.html).

1. How many singly excited configurations do you expect to find there?
2. What is the energy of the first excited state within the single orbital approximation? Is this result in line with the previous exercise?
3. Plot the absorption spectrum for the Cd33Se33 system in the energy interval 0-2 eV using a sigma value of 0.1. (Suggestion: Import the `convolute` function with `from nanoqm.analysis import convolute` and have a look at the script [`convolution.py`](https://github.com/SCM-NV/nano-qmflows/blob/master/scripts/qmflows/convolution.py#L45-L52))

## 3. The distribute_derivative_couplings workflow
<a name="#derivative_couplings"></a> [Back to TOC](#toc)

The last twenty points of a ground state molecular dynamics trajectory for the Cd33Se33 system have been distributed into four chunks, for which the overlaps and couplings have been calculated according to the first two parts of the [Derivative Couplings tutorial](https://qmflows-namd.readthedocs.io/en/latest/derivative_couplings.html#) of nano-qmflows. Follow the tutorial to calculate the overlaps and couplings amongst the missing pairs of points. In your working directory:
- copy the full trajectory `Cd33Se33_MD_last20.xyz`;
- copy the files `chunk_0.hdf5`, `chunk_1.hdf5`, `chunk_2.hdf5`, `chunk_3.hdf5` and merge them into a `chunk_0123.hdf5` file;
- copy locally the file `input.yaml` and customize it with your own path to the merged .hdf5, the full MD trajectory, and the scratch directory;
- create the following `launch.sh` script:

      #!/bin/sh
      #SBATCH --partition=valhalla  --qos=valhalla
      #SBATCH --clusters=faculty
      #SBATCH --account=cyberwksp21
      
      #SBATCH --time=00:10:00
      #SBATCH --nodes=1
      #SBATCH --ntasks-per-node=12
      #SBATCH --mem=32000
       
      eval "$(/projects/academic/cyberwksp21/Software/Conda/Miniconda3/bin/conda shell.bash hook)"
      conda activate qmflows
      module load cp2k/8.1-sse
       
      run_workflow.py -i input.yml
       
and submit your calculation with:
 
    sbatch launch.sh

Use the updated `chunk_0123.hdf5` to retrieve the LUMO-LUMO+1 couplings and plot their value in time.
