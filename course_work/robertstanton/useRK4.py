
import sys
import cmath
import math
import os
from liblibra_core import *
import util.libutil as comn
from libra_py import units
import matplotlib.pyplot as plt
import numpy as np


def main():
    T, Ekin, Epot, Etot, q, p = run_simulations()




def derivatives(Z, params):
    ndof = Z.num_of_rows
    ntraj = Z.num_of_cols
    mass = params["mass"]     # list of dof items
    k = params["force_const"] # -
    q0 = params["q0"]         # -

    q = Z.real()
    p = Z.imag()

    der = CMATRIX(ndof, ntraj)

    for traj in range(ntraj):

        for dof in range(ndof):

            qi =  q.get(dof, traj)
            pi =  p.get(dof, traj)

            val1 = pi/mass[dof]
            val2 = -4*k[dof]*(qi-q0[dof])**3

            dzdt = val1*(1.0+0.0j)+val2*(0.0+1.0j)
            der.set(dof, traj, dzdt )
    return der

def energies(Z, params):
    ndof = Z.num_of_rows
    ntraj = Z.num_of_cols

    mass = params["mass"]     # list of dof items
    k = params["force_const"] # -
    q0 = params["q0"]         # -

    q = Z.real()
    p = Z.imag()

    ekin, epot = 0.0, 0.0

    for traj in range(ntraj):
        for dof in range(ndof):

            qi =  q.get(dof, traj)
            pi =  p.get(dof, traj)

            ekin += 0.5*pi*pi/mass[dof]
            epot += k[dof]*(qi-q0[dof])**4
    ekin = ekin / ntraj
    epot = epot / ntraj

    etot = ekin + epot

    return ekin, epot, etot

def run_simulations():
    
    # Initial conditions
    Z = CMATRIX(1, 1)
    Z.set(0,0, 0.1+0.01j)

    # Potential
    params = {"mass":[1000.0], "force_const":[0.001], "q0":[0.0]}

    # Simulation parameters
    dt = 1.0*units.fs2au
    nsteps = 500

    T, Ekin, Epot, Etot, q, p = [], [], [], [], [], []
    for step in range(nsteps):    
        ekin, epot, etot = energies(Z, params)
    
        T.append(step*dt)
        Ekin.append(ekin)
        Epot.append(epot)
        Etot.append(etot)
        q.append(Z.get(0,0).real)
        p.append(Z.get(0,0).imag)
        #print(F"step = {step} ekin = {ekin} epot = {epot} etot = {etot}")
    
        Z = RK4(Z, dt, derivatives, params)
    print('done')
    return T, Ekin, Epot, Etot, q, p

main()
