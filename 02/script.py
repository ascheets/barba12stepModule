import numpy as np
import matplotlib.pyplot as plt
import sys

#variables to determine discretization 
#in space and time

nx = 41
dx = 2./(nx-1)
nt = int(sys.argv[1])
dt = 0.025

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
        u[i] = un[i] - un[i]*dt/dx*(un[i]-un[i-1])

#plot the results
plt.figure(figsize=(11,7), dpi=100)
plt.plot(np.linspace(0,2,nx),u)
plt.xlim([0,2])
plt.ylim([0,5])
plt.show()
