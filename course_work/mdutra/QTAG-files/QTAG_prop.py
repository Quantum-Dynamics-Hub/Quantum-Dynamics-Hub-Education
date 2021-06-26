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

import QTAG_calc

"""
def deriv_ho():

	for n in range(ntraj):
		qi = qpas.get(n,0)
		pi = #im(grad(psi)/psi)

		val1 = pi/m
		val2 = -k*(qi-x0)

		dzdt = val1*(1.0+0.0j)+val2*(0.0+1.0j)
		deriv.set(0,n,dzdt)

	return(deriv)
"""

def _mom_fit(ntraj,qpas,c):
	mom=MATRIX(ntraj,1);r=MATRIX(ntraj,1)
	gmom=MATRIX(ntraj,1);gr=MATRIX(ntraj,1)
	a=MATRIX(2,2);b=MATRIX(2,2);x=MATRIX(2,2)

	for i in range(ntraj):
		z=complex(0.0,0.0)
		dz=complex(0.0,0.0)
		qpasi=qpas.row(i)
		q1,p1,a1,s1=qpasi.get(0),qpasi.get(1),qpasi.get(2),qpasi.get(3)
		for j in range(ntraj):
			qpasj=qpas.row(j)
			q2,p2,a2,s2=qpasj.get(0),qpasj.get(1),qpasj.get(2),qpasj.get(3)
			term1=np.exp(-0.5*a2*(q1-q2)**2+1.0j*(p2*(q1-q2)+s2))
			z+=c.get(j)*(a2/np.pi)**0.25*term1
			dz-=(a2*(q1-q2)-1.0j*p2)*c.get(j)*(a2/np.pi)**0.25*term1
		ztot=dz/z
		mom.set(i,0,ztot.imag)
		r.set(i,0,ztot.real)

	for m in range(2):
		bb1=0+0j;bb2=0+0j
		for i in range(ntraj):
			x0=qpas.get(i,0)
			z=QTAG_calc.psi(ntraj,qpas,c,x0);zstar=np.conj(z)
			bb1+=mom.get(i)*qpas.get(i,0)**(m)*(z*zstar).real
			bb2+=r.get(i)*qpas.get(i,0)**(m)*(z*zstar).real
		b.set(m,0,bb1.real);b.set(m,1,bb2.real)
		
		for n in range(2):
			aa=0+0j
			for i in range(ntraj):
				x0=qpas.get(i,0)
				z=QTAG_calc.psi(ntraj,qpas,c,x0);zstar=np.conj(z)
				aa+=qpas.get(i,0)**(m+n)*(z*zstar).real
			a.set(m,n,aa.real)

	solve_linsys(a,b,x,1e-1,50000)
	for i in range(ntraj):
		aa=x.get(0,0)+x.get(1,0)*qpas.get(i,0)
		bb=x.get(0,1)+x.get(1,1)*qpas.get(i,0)
		mom.set(i,0,aa);r.set(i,0,bb)
		gmom.set(i,0,x.get(1,0))
		gr.set(i,0,x.get(1,1))
	return(mom,r,gmom,gr)


def symplectic(univ,qpas,c):
	ntraj,dt,mass=univ['ntraj'],univ['dt'],univ['mass']
	qpasn=MATRIX(ntraj,4);an1=MATRIX(ntraj,1)
	mom,r,gmom,gr=_mom_fit(ntraj,qpas,c)

	
	qo=qpas.col(0);ao=qpas.col(2)

	an1.dot_product(ao,gmom)
	an=ao-2.0*an1*dt/mass
#	an=ao
	qn=qo+dt*mom/mass

	for i in range(ntraj):
		qpasn.set(i,0,qn.get(i))
		qpasn.set(i,1,mom.get(i))
		qpasn.set(i,2,an.get(i))
		qpasn.set(i,3,qpas.get(i,3))

	return(qpasn)

