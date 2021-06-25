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
import QTAG_calc

def hamiltonian(ntraj,qpas,ov,nsurf):
	H = CMATRIX(ntraj,ntraj)
	q1=MATRIX(1,1); q2=MATRIX(1,1)
	p1=MATRIX(1,1); p2=MATRIX(1,1)
	a1=MATRIX(1,1); a2=MATRIX(1,1)
	s1=MATRIX(1,1); s2=MATRIX(1,1)
	
	for i in range(ntraj):
		for j in range(ntraj):
			qpasi=qpas.row(i)
			qpasj=qpas.row(j)
			q1.set(0,0,qpasi.get(0));p1.set(0,0,qpasi.get(1));a1.set(0,0,qpasi.get(2));s1.set(0,0,qpasi.get(3))
			q2.set(0,0,qpasj.get(0));p2.set(0,0,qpasj.get(1));a2.set(0,0,qpasj.get(2));s2.set(0,0,qpasj.get(3))
			ke=(-1.0**2/(2.0*mass))*gwp_kinetic(q1,p1,s1,a1/2,q2,p2,s2,a2/2)

			v=QTAG_pots.LHA(qpasi,qpasj,nsurf)
			H.set(i,j,ke+v*ov.get(i,j))
	return(H)

def coupling(ntraj,qpas1,qpas2,ov12,nsurf):
	H = CMATRIX(ntraj,ntraj)

	for i in range(ntraj):
		for j in range(ntraj):
			qpasi=qpas1.row(i)
			qpasj=qpas2.row(j)

			v=QTAG_pots.LHA(qpasi,qpasj,nsurf)
#			v2=QTAG_pots.exact_gauss_cpl(qpasi,qpasj)
#			print(i,j,np.abs(v*ov12.get(i,j)),np.abs(v2*ov12.get(i,j)),np.abs(v-v2),sep=' ',end='\n',file=f_vtest)
			H.set(i,j,v*ov12.get(i,j))
	return(H)

def basis_diag(m,H,ov,b):
	evals=CMATRIX(m,m)
	evecs=CMATRIX(m,m)
	solve_eigen(H,ov,evals,evecs,0)
	ct=evecs.T().conj()
	c_new=evecs*(exp_(evals,-dt*1.0j))*ct*b

	return(c_new)

def _nonad_assemble(chk,m,n,matrix1,matrix2,matrix3):
	if chk=='real':
		if m==n:
			mtot=MATRIX(2*m,2*n)
		else:
			mtot=MATRIX(2*m,n)
	elif chk=='cplx':
		if m==n:
			mtot=CMATRIX(2*m,2*n)
		else:
			mtot=CMATRIX(2*m,n)

	listm=[]; listn=[]
	for i in range(m):
		listm.append(i)
	for j in range(n):
		listn.append(j)

	if n==m:
		listm=[]; listn=[]
		for i in range(m):
			listm.append(i)
		for j in range(n):
			listn.append(j)

		push_submatrix(mtot,matrix1,listm,listn)
		push_submatrix(mtot,matrix2,[m+lm for lm in listm],[n+ln for ln in listn])
		push_submatrix(mtot,matrix3,listm,[n+ln for ln in listn])
		push_submatrix(mtot,matrix3,[m+lm for lm in listm],listn)

	else:
		listm=[]; listn=[]
		for i in range(m):
			listm.append(i)
		for j in range(n):
			listn.append(j)

		push_submatrix(mtot,matrix1,listm,listn)
		push_submatrix(mtot,matrix2,[m+lm for lm in listm],listn)

	return(mtot)

f_traj1=open('gbc1.txt', 'w')
f_traj2=open('gbc2.txt', 'w')
f_obs=open('obs.txt', 'w')
f_vtest=open('vtest.txt','w')

univ, wf0=QTAG_init.params()
ntraj,a0,rcut=univ['ntraj'],univ['a0'],univ['rcut']
niter,dt,mass=univ['niter'],univ['dt'],univ['mass']
nout,t12_chk=univ['nout'],univ['t12_chk']
qpas1=MATRIX(ntraj,4);qpas2=MATRIX(ntraj,4)
b1=CMATRIX(ntraj,1);b2=CMATRIX(ntraj,1)
c1_new=CMATRIX(ntraj,1);c2_new=CMATRIX(ntraj,1)
ct_old=CMATRIX(2*ntraj,1);ct_new=CMATRIX(2*ntraj,1)

dummy=CMATRIX(ntraj,ntraj)

pop_list=[]
for i in range(ntraj):
	pop_list.append(i)

qpas1=QTAG_init.grid(univ,wf0,qpas1)
qpas2=QTAG_init.grid(univ,wf0,qpas2)

b1=QTAG_init.coeffs(ntraj,wf0,qpas1,1)
b2=QTAG_init.coeffs(ntraj,wf0,qpas2,2)

bt=_nonad_assemble("cplx",ntraj,1,b1,b2,dummy)
#qpast=_nonad_assemble("real",ntraj,4,qpas1,qpas2,dummy)

t=0.0

for iter in range(niter):
#	if iter%nout==0:
#		print(iter)

	ov1=QTAG_calc.overlap(ntraj,qpas1,qpas1)
	ov2=QTAG_calc.overlap(ntraj,qpas2,qpas2)
	ov12=QTAG_calc.overlap(ntraj,qpas1,qpas2)

	H1=hamiltonian(ntraj,qpas1,ov1,1)
	H2=hamiltonian(ntraj,qpas2,ov2,2)
	Hcpl=coupling(ntraj,qpas1,qpas2,ov12,3)

	ovt=_nonad_assemble("cplx",ntraj,ntraj,ov1,ov2,dummy)
	Ht=_nonad_assemble("cplx",ntraj,ntraj,H1,H2,Hcpl)

	ct_new=basis_diag(2*ntraj,Ht,ovt,bt)
	pop_submatrix(ct_new,c1_new,pop_list,[0])
	pop_submatrix(ct_new,c2_new,[ntraj+i for i in pop_list],[0])

	qpas1n=QTAG_prop.symplectic(univ,qpas1,c1_new)
	ov_no=QTAG_calc.overlap(ntraj,qpas1n,qpas1)
	b1=ov_no*c1_new


	norm2=QTAG_calc.norm(c2_new,ov2).real
	if norm2 < t12_chk:
#		for i in range(ntraj):
#			for j in range(4):
#				qpas2n.set(i,j,qpas1n.get(i,j))
		qpas2n=qpas1n
	else:
		qpas2n=QTAG_prop.symplectic(univ,qpas2,c2_new)

	ov_no=QTAG_calc.overlap(ntraj,qpas2n,qpas2)
	b2=ov_no*c2_new

	qpas1=QTAG_calc.update(ntraj,4,qpas1,qpas1n)
	qpas2=QTAG_calc.update(ntraj,4,qpas2,qpas2n)

	bt=_nonad_assemble("cplx",ntraj,1,b1,b2,dummy)
	ct_old=QTAG_calc.update(ntraj,1,ct_old,ct_new)

	if iter%nout==0:
		out1=list(qpas1.get(i,0) for i in range(ntraj))
		out2=list(qpas2.get(i,0) for i in range(ntraj))
		print(t,*out1,sep=' ',end='\n',file=f_traj1)
		print(t,*out2,sep=' ',end='\n',file=f_traj2)
		print(t,QTAG_calc.norm(c1_new,ov1).real, QTAG_calc.norm(c2_new,ov2).real, \
		      QTAG_calc.energy(ct_new,Ht).real,sep=' ', end='\n', file=f_obs)
	t+=dt

QTAG_calc.wf_print(ntraj,qpas1,c1_new,'wft1.txt')
QTAG_calc.wf_print(ntraj,qpas2,c2_new,'wft2.txt')
