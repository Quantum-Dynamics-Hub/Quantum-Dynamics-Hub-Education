import sys
import cmath
import math
import os
from liblibra_core import *
import util.libutil as comn

#import libra_py
from libra_py import data_outs
import matplotlib.pyplot as plt
import numpy as np


def params():
	global wf0; global univ
	wf0 = {"q":-2.0, "p":0.0, "a":1.0, "s":0.0}
	univ = {"ntraj":25, "dt":0.01, "rcut":1e-12, "a0":10.0, "niter":800, "mass":1.0, "t12_chk":0.9, "nout":1}

	return(univ,wf0)

def grid(univ,wf0,qpas):
	ntraj=univ['ntraj']
	rcut=univ['rcut']
	a0=univ['a0']

	xlow=wf0['q']-np.sqrt(-0.5/wf0['a']*np.log(rcut))
	xhi=wf0['q']+np.sqrt(-0.5/wf0['a']*np.log(rcut))
	qs=np.linspace(xlow,xhi,num=ntraj)

	for i in range(ntraj):
		qpas.set(i,0,qs[i])
		qpas.set(i,1,wf0['p'])
		qpas.set(i,2,wf0['a']*a0)
		qpas.set(i,3,0.0)

	return(qpas)


def coeffs(ntraj,wf0,qpas,nsurf):
	b=CMATRIX(ntraj,1)

	if nsurf==1:
		for i in range(ntraj):
			qpasi=qpas.row(i)
			q1,p1,a1,s1=qpasi.get(0),qpasi.get(1),qpasi.get(2),qpasi.get(3)
			q2,p2,a2,s2=wf0['q'],wf0['p'],wf0['a'],wf0['s']
		
			b.set(i,0,gwp_overlap(q1,p1,s1,a1/2,q2,p2,s2,a2))
	elif nsurf==2:
		for i in range(ntraj):
			b.set(i,0,0+0j)
	return(b)

