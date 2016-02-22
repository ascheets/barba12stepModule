import numpy as np
import pylab
import matplotlib.pyplot as plt
import time, sys

#returns an error when running from terminal
#%matplotlib inline

#number of spatial steps
nx = 41
#size of each spatial step
dx = 2./(nx-1)
#number of time steps
nt = 25
#length of each time step
dt = 0.025
#wavespeed
c = 1

#initialize u -> velocity
u = np.ones(nx)
#initial condition -> square wave
u[0.5/dx : 1/dx+1] = 2
print(u)

plt.plot(np.linspace(0,2,nx), u)
pylab.show()

#initialize temp array
#array will be used for storing u at time "n"
un = np.ones(nx)

#nested for loop iterates through all points in space
#for each moment in time

#iterate through all time steps
for n in range(nt):
    #copy old version of u
    un = u.copy()
    #range(stop) defaults to zero
    #iterate over spatial domain
    for i in range(1,nx):
        #expression relates to forward difference time, 
        #backward difference space
        u[i] = un[i] - c*dt/dx*(un[i]-un[i-1])

plt.plot(np.linspace(0,2,nx),u)
pylab.show()
