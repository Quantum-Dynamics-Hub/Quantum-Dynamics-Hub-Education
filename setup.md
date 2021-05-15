---
title: Setup
---

> ## Warning
>
> This webpage is still under development. Please check for updates.
>
{:.discussion}


## Tools you will need
1. Getting to the UB network (account set up)
2. Open On Demand
3. SEAGrid
4. Libra installation
5. Text Editors


## 1. Setting up your UB CCR account and getting access to the UB CCR supercomputer and infrastructure
Below are the steps you’ll need to take in order to connect to our resources:

### Step 1.

In order to connect to our machines, you have to be on the UB network.

You'll need to download and install UB's Cisco VPN client. You can download the client software from here:
[For Windows](http://www.buffalo.edu/ubit/service-guides/software/downloading/windows-software/managing-your-software/anyconnect.html), 
[For MacOS](http://www.buffalo.edu/ubit/service-guides/software/downloading/macintosh-software/managing-mac-software/anyconnect.html) 

You **MUST** use the Cisco AnyConnect client with this account, please ignore any information about the FortiClient VPN software.

*The login information to download and connect to the UB VPN client will be given to you over emain, when your are accepted*

NOTE: When you start the Cisco software the first time you will need to enter the following in the "Connect to:" box: vpn.buffalo.edu
and then select UBVPN from the group drop down menu.

### Step 2 (one time only).

You will receive a separate email containing a link to our user portal to setup your account and change your password. 

You will need to already be connected to the UB VPN to get this site and to access the one-time link. 

The link contained in this email only lasts for 24 hours so if you're not able to get this done before it expires, you can go to:
  https://idm.ccr.buffalo.edu and click the "forgot your password?" link to generate a new one. 
  Your CCR username is the same as the one you’re using for the VPN (given to you in the email), but these accounts are not connected in any other way.
 


## 3. Terminal

The terminal is an interface in which you can type and execute text based commands. It is important to use the terminal to 
run many computational chemistry software packages. There are several different types of terminal interfaces, called shells.
In this tutorial, we will focus on using one of the most common shells, the bash shell. 
How you acquire a bash shell terminal depends on the type of computer you have.

### Linux
If you are using a Linux computer, you probably already know how to open the terminal window. 
If the Terminal is not shown in menu of programs, you can use the key combination CTRL + ATL + T to open the terminal window.

### Mac OS X
On Mac OS X, a Terimanl application is built into your system. Open the Terminal from Applications -> Utilities -> Terminal.

### Windows
Windows has a built in command line interface. To access it, click the Windows Key + R, type cmd, press Enter.
**However,** this interface is not a bash application. Therefore, the commands for navigating and creating files discussed below
 will not be the same. We recommend you installing the [Windows Subsystem for Linux]()


### Resources
- [Using the command line](https://ryanstutorials.net/linuxtutorial/commandline.php)
- [Navigation in bash](https://ryanstutorials.net/linuxtutorial/navigation.php)
- [Making and removing directories, copying and deleting files](https://ryanstutorials.net/linuxtutorial/filemanipulation.php)


## 3. Python & Jupyter & Matplotlib

All these tools are already available via VIDIA, but if you prefer to work on your local machine, you may simply install the Anaconda python.

[Python](https://python.org/) is a popular language for scientific computing, and great for general-purpose programming as well. 
Installing all of its scientific packages individually can be a bit difficult, however, so we recommend the all-in-one installer Anaconda.
Regardless of how you choose to install it, *please make sure you install Python version 3.x (e.g., 3.4 is fine, 2.7 is not)*.

If you aren't very comfortable with Python yet, [this resource](https://www.tutorialspoint.com/python/index.htm) could be a good starting point.

[Jupyter](https://jupyter.org/) is an tool to enble interactive experience with Python and any other packages that can be called via Python

[Matplotlib](https://matplotlib.org/) is a plotting library than can be called by Python. When integrated into Jupyter notebooks, it allows you
to plot your results on the go - as soon as you obtain it. Please check the official Matplotlib site for a 
[great collection of examples](https://matplotlib.org/gallery/index.html)

Here is a link of more specific topics on Jupyter
* How to run Jupyter: [link](https://nbviewer.jupyter.org/github/jupyter/notebook/blob/master/docs/source/examples/Notebook/Running%20Code.ipynb)
* How to plot 2D and 3D figures: [link](https://nbviewer.jupyter.org/github/jrjohansson/scientific-python-lectures/blob/master/Lecture-4-Matplotlib.ipynb)
* A gallery of interesting Jupyter examples: [link](https://github.com/jupyter/jupyter/wiki/A-gallery-of-interesting-Jupyter-Notebooks#introductory-tutorials)


## 4. Text Editors

You will often need to create or read text files.  Opening a text file in a word processing program,
like Microsoft Word or Google Docs, introduces a lot of formatting that is not needed. 
You need to use a text editor to read and write these files. There are many choices. You don't need to learn to 
use all of these at the beginning, just find one that works for you. 

### Vi/vim
Vi/vim is one of the most ubiquitous text editors. It is installed on virtually every Linux computer in the world, 
so if you ever log on to a unfamiliar machine, it will be available to you. 
Vi is accessed from the command line; it doesn't have or need a graphical interface so it can operate on the most bare bones computers.
However, it is not intuitive to use and can be difficult for beginners.
- [Interactive vim tutorial](https://www.openvim.com/)
- [Getting started with vi](https://ryanstutorials.net/linuxtutorial/vi.php)

### Far 3
[Far3](https://www.farmanager.com/) is a very handy file and archive manager. Features rather convenient graphical interface and some
convenient manipulation options such as cutting any blocks of text, syntax highlightling, remote drive access, etc. Don't get
scared by its old-style "Norton/Midnight Commander" look. You'll love it. If you run Windows, this is definitely our recommendation. 

### Atom
[Atom](https://atom.io/) is a modern text editor that is very intuitive to use. You probably don't even need 
to read the tutorial below to figure out how to create and save files. Standard downloads are available for Linux, Mac, and Windows.
- [Getting started with atom](https://flight-manual.atom.io/getting-started/sections/atom-basics/)

### Emacs
Similar to vi/vim, emacs is a command line text editor that is already part of almost all Linux distributions.
Emacs can also be used as an RSS reader or file manager.
- [Emacs tutorial](https://www.gnu.org/software/emacs/tour/) The section called Beginner tips near the bottom of the page is particularly helpful.

### Sublime
[Sublime](https://www.sublimetext.com/) is an intuitive text editor that makes looking at files with multiple sections easy. 
t allows split screen editing and is very customizable.



---

> ## Prerequisites
> * the basic experience with Python language
> * the basic knowledge of terminal
>
{:.prereq}




{% include links.md %}
