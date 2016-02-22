#Convergence and the CFL Condition
import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
import sys

def linearconv(nx):
    dx = 2./(nx-1)
    nt = 20 #int(sys.argv[1])
    dt = 0.025
    c = 1
    
    #initialize array
    u = np.ones(nx)
    #set u = 2 between 0.5 and 1
    u[0.5/dx : 1/dx+1] = 2

    #initialize placeholder array un
    un = np.ones(nx)

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

#plot the results
plt.plot(np.linspace(0,2,nx),u)
pl.show()
