import matplotlib.pyplot as plt
import numpy as np
from math import *

N_integral = 120
dt = 0.05
omega = 2*3.14*1.0
alpha = 1.0

def filter(k):
    return exp(-alpha*k*dt) * (cos(omega*k*dt)) 

N0=100
x = np.linspace(0,1,4*N0)
y = np.concatenate([np.zeros(3*N0),np.ones(N0),np.ones(N0),4*np.ones(N0)])

stored = []

normalization = 0
for k in range(N_integral):
    normalization = normalization+filter(k)

z = np.zeros(len(y))
for k_z in range(len(y)-N_integral):
    for k in range(N_integral):
        z[len(y)-1-k_z] = z[len(y)-1-k_z] + y[len(y)-1-k_z-k] * filter(k) #* (+cos(omega*k*dt)-alpha/omega*sin(omega*k*dt)) * 2# / N_integral


plt.plot(y,'r')
plt.plot(y+(z/normalization-y)*0.2,'b')
#plt.plot(stored,'b')
plt.show()

