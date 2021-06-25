import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

class QT(object):
    def __init__(self, x, S, S1, S2, C, C1, C2, a, t=0.):
        self.x, self.S, self.S1, self.S2, self.C, self.C1, self.C2, self.a = map(np.array, (x, S, S1, S2, C, C1, C2, a))
        self.xdot = np.zeros(self.x.shape)
        self.t = t

    # velocity of outer points
    def dr(self, r, t):
        return self.S1[::len(self.x)-1]

    # C-amplitude time-derivative
    def dC(self, c, t):
        (c1, c2) = self.mls_d(self.xo, c)
        return (self.xdot-self.S1)*c1 - self.S2/2.

    # phase time-derivative
    def dS(self, s, t):
        (s1, s2) = self.mls_d(self.xo, s)
        return self.xdot*s1 - .5*s1**2 + .5*(self.C1**2 + self.C2) - self.v()

    # classical potential
    def v(self):
        return np.poly1d(a)(self.xo)

    # derivative calculator
    def mls_d(self,x,y):
        k=4
        N=len(x)
        M=k+1
        dyf, ddyf, E, W, B = map(np.zeros, (N, N, N, (N,N), (N,M)) )
        a = 1.
        for k in range(N):
            for i in range(N):
                W[i][i] = np.exp(-a*(x[i]-x[k])**2)
                for j in range(M):
                    B[i][j] = (1./math.factorial(j))*(x[i]-x[k])**j
            BtBi = np.linalg.pinv( B.T @ (W @ B) )
            c = BtBi @ (B.T @ (W @ y))
            dyf[k] = c[1]
            ddyf[k] = c[2]
        return (dyf,ddyf)


    # time stepping
    def time_step(self, dt, tspan, N):
        self.solx = np.zeros((len(tspan),len(self.x))) #holds solution for plotting
        for i in range(len(tspan)):
            self.xo = self.x # holds current positions
            self.x[::len(self.x)-1] = odeint(self.dr, self.x[::len(self.x)-1], [0,dt])[-1] # update outer points
            self.x = np.linspace(self.x[0], self.x[-1], self.x.size) # update inner points
            self.xdot[::len(self.x)-1] = self.dr(self.x[::len(self.x)-1], dt) # outer grid velocities
            self.xdot[1:-1] = (self.x[1:-1] - self.xo[1:-1])*(1./dt) # inner grid velocities
            self.S = odeint(self.dS, self.S, [0, dt])[-1] # update phase
            self.C = odeint(self.dC, self.C, [0, dt])[-1] # update C-amplitude
            (self.C1, self.C2) = self.mls_d(self.x, self.C) # obtain spatial-derivatives of C-amplitude
            (self.S1, self.S2) = self.mls_d(self.x, self.S) # spatial derivatives of phase
            self.solx[i,:] = self.x
            self.t = tspan[i]
            print('Time : ', tspan[i])
            print('Norm :', (self.x[-1]-self.x[0])/N * np.sum(np.exp(2*self.C)))

################################################################

# gaussian wave packet
def gauss_x(b, x, x0):
    return (2. * b / np.pi)**.25 * np.exp(-b * (x - x0)**2)

################################################################


# potential parameters
a = np.array([0., 0., 0.])

b0 = .5 # width
x0 = .0 # displacement
wn = 3.0 # wavenumber

# amplitude
N = 25 # number of fluid elements

# initial conditions
x = np.linspace(x0-2, x0+2, N) # support
R = gauss_x(b0, x, x0) # amplitude

# C-amplitude
C = np.log(R)
C2 = np.zeros(x.shape)-2.*b0
C1 = C2*(x-x0)

# phase
S2 = np.zeros(x.shape)
S1 = S2 + wn
S = S1 * x

Q = QT(x=x,
        S=S,
        S1=S1,
        S2=S2,
        C=C,
        C1=C1,
        C2=C2,
        a=a)

#time
dt = .005
tspan = np.arange(0.0, 5.0, dt)

# run time-steps for free packet
Q.time_step(dt,tspan,N)

plt.plot(tspan, Q.solx)
plt.ylabel('x')
plt.xlabel('t')
plt.show()

# potential parameters
a = np.array([2., 0., 0.])

# run time-steps for harmonic oscillator
Q.time_step(dt,tspan,N)

plt.plot(tspan, Q.solx)
plt.ylabel('x')
plt.xlabel('t')
plt.show()
