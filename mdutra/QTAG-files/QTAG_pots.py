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

def model_ho(x,nsurf):
	k1,k2=10.0,10.0
	x0=1.0;y0=5.0*np.sqrt(k1)
	d1=1.0;d2=0.5*np.sqrt(k1);d3=2.0
	if (nsurf==1):
		fx=0.5*k1*x**2
		dfx=k1*x
		d2fx=k1
	elif (nsurf==2):
		fx=0.5*k2*(x-x0)**2+y0
		dfx=k2*(x-x0)
		d2fx=k2
	elif (nsurf==3):
		fx=d1*np.exp(-d2*(x-d3)**2)
		dfx=-2.0*d1*d2*(x-d3)*np.exp(-d2*(x-d3)**2)
		d2fx=-2.0*d1*d2*np.exp(-d2*(x-d3)**2)-dfx*2.0*d2*(x-d3)
	return(fx,dfx,d2fx)

def model_t1(x,nsurf):
	a=0.01;b=1.147
	d1=0.005;d2=1.0;d3=0.0
	if (nsurf==1):
		fx=a*(1.0+np.tanh(b*x))
		dfx=a*b*(1.0-np.tanh(b*x)**2)
		d2fx=-2.0*a*b**2*np.tanh(b*x)*(1.0-np.tanh(b*x)**2)
	elif (nsurf==2):
		fx=a*(1.0-np.tanh(b*x))
		dfx=-a*b*(1.0-np.tanh(b*x)**2)
		d2fx=2.0*a*b**2*np.tanh(b*x)*(1.0-np.tanh(b*x)**2)
	elif (nsurf==3):
		fx=d1*np.exp(-d2*(x-d3)**2)
		dfx=-2.0*d1*d2*(x-d3)*np.exp(-d2*(x-d3)**2)
		d2fx=-2.0*d1*d2*np.exp(-d2*(x-d3)**2)-dfx*2.0*d2*(x-d3)

	return(fx,dfx,d2fx)

def exact_gauss_cpl(qpasi,qpasj):
	v=complex(0.0,0.0)
	k1,k2=10.0,10.0
	d1=1.0;d2=0.5*np.sqrt(k1);d3=2.0
#	d1=0.005;d2=1.0;d3=0.0
	q1,p1,a1=qpasi.get(0), qpasi.get(1), qpasi.get(2)
	q2,p2,a2=qpasj.get(0), qpasj.get(1), qpasj.get(2)

	et1=-(q1-d3)**2*a1**2/2.0-(q2-d3)**2*a2**2/2.0+(p1-p2)**2/2.0
	et2=(q1-d3)*((d3-q2)*a2+1j*(p1-p2))*a1
	et3=1j*(q2-d3)*(p1-p2)*a2
	et4=(a1+2*d2+a2)*(a1+a2)

	v=np.exp(2.0*d2*(et1+et2+et3)/et4)*d1*np.sqrt(2.0*a1+2.0*a2)/np.sqrt(2.0*a1+4.0*d2+2.0*a2)

	return(v)

def LHA(qpasi,qpasj,nsurf):
	v=complex(0.0,0.0)
	q1,p1,a1=qpasi.get(0), qpasi.get(1), qpasi.get(2)
	q2,p2,a2=qpasj.get(0), qpasj.get(1), qpasj.get(2)

	z=complex(a1*q1+a2*q2,p2-p1)/(a1+a2)
	xs=[q1,q2]
	for x in xs:
		vx,dvx,d2vx=model_ho(x,nsurf)
		vv0=vx-dvx*x+d2vx/2.0*x**2
		vv1=-d2vx*x+dvx
		vv2=d2vx/2.0
		v+=0.5*(vv0+vv1*z+vv2*(z**2+1.0/(a1+a2)))

	return(v)


