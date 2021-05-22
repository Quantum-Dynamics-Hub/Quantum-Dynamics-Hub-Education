---
title: "Introduction to our CyberInfrastructure. Setups and key links."
date: June 14, 2021, 11:00 am - 5:00 pm EDT
---

## Tools you will need
1. Overview of the UB CyberInfrastructure
2. [UB OnDemand](https://ubccr.freshdesk.com/support/solutions/articles/13000039875-ccr-ondemand-portal)
3. Using Jupyter on the OnDemand gateway

### 1. Overview of the UB CyberInfrastructure

In addition to your home directory, you have access to the workshop directory:

    /projects/academic/cyberwksp21

This folder contains the following sub-directories:

* `Instructors_material`
    Contains the examples or working data for the tutorials. You can access it to copy the content you need, but
    **Do not edit or view files in this directory at any time**. You can still `ls` directories to see the content, but
    do not `vi` the files. You may accidentally edit or change the files, which may affect other users, so please be mindful. 

* `Modules`
    Contains definitions and setups of the environmental variables for all users. **Do not edit or view files in this directory at any time**
    If you are curious on what is in there, feel free to first copy it to your working directory and then explore.

* `Software`
    Contains the installations of some packages, such as Conda (with the corresponding environments), or Libra.
    **Do not edit or view files in this directory at any time**. You can `ls` and `cd` there, to explore the content of
    directories, but do not operate on files. If you are curious about some of them, make your own copy first. 
    
* `Students`
    First, please create your own directory in this folder.

    This will be your working directory (apart from your home directory). This is where you can keep your data 
    and run some (small) calculations. Data in this directory are shared among the participants (those who have
    access to the UB resources), so you can take the advantage of that: e.g. if something doesn't work, you may 
    check our outher students' directories, but again - **Do not edit or view files in other students' directories at any time**


          
### 2. Introduction to the OnDemand
A very detailed introduction is given [in this video](https://ub.hosted.panopto.com/Panopto/Pages/Viewer.aspx?id=c5c088f6-ba8c-4210-8d87-ab9f0104f54e)

[This website](https://ubccr.freshdesk.com/support/solutions/articles/13000039875-ccr-ondemand-portal) also 
provides the ddetailed step-by-step instructions for logging into the system, as well as a general overview of
the available tools.

Before you can use our Python installations via Jupyter, you need to edit your `.bashrc` file in your home directory.

1. Once you have logged into OnDemand, go to Clusters -> Faculty Cluster Shell.

2. Go to home directory and edit the `.bashrc` file

        cd 
        vi .bashrc

3. Add the following line in the file

        module use /projects/academic/cyberwksp21/Modules

This will enable both Jupyter (called from the OnDemand) access the 
right installations as well as will enable you to acccess some specialized 
modules. 

4. Make sure you restart the terminal or `source .bashrc` for the above change to take effect

5. To check which modules are available:

        module avail 

As the result, you shall be able to see modules that are installed system-wide as well as specific modules for the
workshop (defined in `/projects/academic/cyberwksp21/Modules` )


![](../fig/1_episode/modules.png){:width="80%"}

### 3. Using Jupyter on the OnDemand gateway

1. Make sure you include the above `module use /projects/academic/cyberwksp21/Modules` line in
   your `.bashrc` file

2. When you start Jupyter notebook you will have access to only your home directory.
   It is advisable to keep the working file in the /projects/academic/cyberwksp21/Students/\<yourid\>

   To acceess such a directory, create a soft link in your home directory:

    ln -s /projects/academic/cyberwksp21 ~/workshop

3. When you create new or open an existing Jupyter notebook, make sure to select correct kernel:

        Kernels -> Change kernel -> Python 3 (libra-latest)

 

{% include links.md %}

