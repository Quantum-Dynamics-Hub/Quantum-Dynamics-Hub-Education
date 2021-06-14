---
title: "Coding Molecular Dynamics. Intro to Libra."
date: June 14, 2021, 3:40 pm - 5:00 pm EDT

---

<a name="toc"></a>
# Table of Content
1. [Molecular dynamics](#md)
2. [Introduction to Libra](#intro_libra)


<a name="md"></a>
## 1. Molecular Dynamics
[Back to TOC](#toc)

### 1.1. [Presentation](/files/episode_1/Akimov-MD.pdf)

### 1.2. Coding Exercises 

  * Euler MD
  * Leap-frog (coordinate Verlet)
  * Velocity verlet

### 1.3. Libra exmaples

  * [Runge-Kutta for classical MD](https://github.com/compchem-cybertraining/Tutorials_Libra/tree/master/2_integrators/1_runge_kutta_4th_order)

  * Lattice MD

     `cd /projects/academic/cyberwksp21/Instructors_material/alexeyak/libra/demos_from_libra_code/notebooks/Example4_md/Example1_lattice`

  * Atomistic MD

     `cd /projects/academic/cyberwksp21/Instructors_material/alexeyak/libra/demos_from_libra_code/tests/test_md`



<a name="intro_libra"></a>
## 2. Introduction to Libra
[Back to TOC](#toc)

### 2.1. [Presentation](/files/episode_1/Akimov-Libra-overview.pdf)

### 2.2. Get your copy of tutorials

    git clone https://github.com/compchem-cybertraining/Tutorials_Libra.git

   or 

    cd /projects/academic/cyberwksp21/Instructors_material/alexeyak/libra


### 2.3. Examples 

   * Discuss data types

   * [Runge-Kutta for quantum dynamics](https://github.com/compchem-cybertraining/Tutorials_Libra/tree/master/2_integrators/2_runge_kutta_4_for_Liouville)

   * Molecular integrals
 
         cd /projects/academic/cyberwksp21/Instructors_material/alexeyak/libra/demos_from_libra_code/tests/test_libint2

   * Built-in Electronic structure calculations: 

         cd /projects/academic/cyberwksp21/Instructors_material/alexeyak/libra/Tutorials_Libra/5_electronic_structure/1_eht/1_compact1

   * Kinetic MC example to model diffusion

         cd /projects/academic/cyberwksp21/Instructors_material/alexeyak/libra/demos_from_libra_code/notebooks/Example5_diffusion

   * random number sampling - section 7.4

   * crystal builder - section 12.1 (may be useful for the nano-qmflows/CAT tutorials)

   * Molecular rotation

         cd /projects/academic/cyberwksp21/Instructors_material/alexeyak/libra/demos_from_libra_code/tests/example_8_rotation

     Exercise: can you combine examples showing py3Dmol, imageio, and the current rotation script to produce the animateg gifs of the
     rotation?









