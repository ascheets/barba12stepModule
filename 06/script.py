#array operations with numpy
import numpy as np

u = np.array((0,1,2,3,4,5))

for i in range(1,len(u)):
    print u[i]-u[i-1]

u[1:]-u[0:-1]

