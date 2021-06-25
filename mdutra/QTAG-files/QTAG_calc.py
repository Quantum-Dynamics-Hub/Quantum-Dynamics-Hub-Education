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

import QTAG_init
import QTAG_pots
import QTAG_prop

def psi(ntraj,qpas,c,x0):
	wf=0+0j
	for i in range(ntraj):
		q1,p1,a1,s1=qpas.get(i,0),qpas.get(i,1),qpas.get(i,2), qpas.get(i,3)
		wf+=c.get(i)*(a1/np.pi)**0.25*np.exp(-a1/2.0*(x0-q1)**2+1j*(p1*(x0-q1)+s1))
	return(wf)

def energy(c,H):
	e=c.T().conj()*H*c
	return(e.get(0))

def norm(c,ov):
	n=c.T().conj()*ov*c
	return(n.get(0))

def update(ntraj,ndim,a,b):
	for n in range(ndim):
		for i in range(ntraj):
			a.set(i,n,b.get(i,n))
	return(a)

def overlap(ntraj,qpas1,qpas2):
	ov_mat=CMATRIX(ntraj,ntraj)

	for i in range(ntraj):
		for j in range(ntraj):
			qpasi=qpas1.row(i)
			qpasj=qpas2.row(j)
			q1,p1,a1,s1=qpasi.get(0),qpasi.get(1),qpasi.get(2), qpasi.get(3)
			q2,p2,a2,s2=qpasj.get(0),qpasj.get(1),qpasj.get(2), qpasj.get(3)

			ov_mat.set(i,j,gwp_overlap(q1,p1,s1,a1/2,q2,p2,s2,a2/2))
	return(ov_mat)

def wf_print(ntraj,qpas,c,filename):
	f_wf=open(filename,'w')

	for x0 in np.linspace(-8.0,8.0,num=1000):
		z=psi(ntraj,qpas,c,x0);zstar=np.conj(z)
		print(x0,np.abs(z),(zstar*z).real,sep=' ', end='\n', file=f_wf)

	f_wf.close()
	return()

