from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

#variable declarations
nx = 41
ny = 41
nt = 480
c = 1.2
dx = 2.0/(nx-1)
dy = 2.0/(ny-1)
sigma = 0.0009
nu = 0.01
dt = sigma*dx*dy/nu

x = np.linspace(0,2,nx)
y = np.linspace(0,2,ny)

#num y corresponds to rows
#num x corresponds to cols
u = np.ones((ny,nx))
v = np.ones((ny,nx))
un = np.ones((ny,nx))
vn = np.ones((ny,nx))
comb = np.ones((ny,nx))

#assign initial conditions
#u[y region, x region]
u[0.5/dy:1./dy+1., 0.5/dx:1./dx+1.] = 10.
v[0.5/dy:1./dy+1., 0.5/dx:1./dx+1.] = 1.

#plot ICs
fig = plt.figure(figsize=(11,7), dpi=100)
ax = fig.gca(projection='3d')
X,Y = np.meshgrid(x,y)
wire1 = ax.plot_wireframe(X,Y,u[:], cmap=cm.coolwarm)
wire2 = ax.plot_wireframe(X,Y,v[:], cmap=cm.coolwarm)
plt.show()

#advance in time
for n in range(nt+1):
    #store copy of old values
    un = u.copy()
    vn = v.copy()

    #code the numerical scheme, array operations
    u[1:-1,1:-1] = un[1:-1,1:-1] - dt/dx*un[1:-1,1:-1]*(un[1:-1,1:-1]-un[0:-2,1:-1]) \
                   - dt/dy*vn[1:-1,1:-1]*(un[1:-1,1:-1]-un[1:-1,0:-2]) \
                   + nu*dt/dx**2*(un[2:,1:-1] - 2.*un[1:-1,1:-1] + un[0:-2,1:-1]) \
                   + nu*dt/dy**2*(un[1:-1,2:] - 2.*un[1:-1,1:-1] + un[1:-1,0:-2])

    v[1:-1,1:-1] = vn[1:-1,1:-1] - dt/dx*un[1:-1,1:-1]*(vn[1:-1,1:-1]-vn[0:-2,1:-1]) \
                   - (dt/dy)*vn[1:-1,1:-1]*(vn[1:-1,1:-1]-vn[1:-1,0:-2]) \
                   + nu*dt/dx**2*(vn[2:,1:-1] - 2.*vn[1:-1,1:-1] + vn[0:-2,1:-1]) \
                   + nu*dt/dy**2*(vn[1:-1,2:] - 2.*vn[1:-1,1:-1] + vn[1:-1,0:-2])

    #reestablishing bcs
    u[0,:] = 1 #first line across x
    u[-1,:] = 1 #last line across x
    u[:,0] = 1 #first line down y
    u[:,-1] = 1 #last line down y

    v[0,:] = 1
    v[-1,:] = 1
    v[:,0] = 1
    v[:,-1] = 1

#plotting
fig = plt.figure(figsize=(11,7), dpi=100)
ax = fig.gca(projection='3d')
X,Y = np.meshgrid(x,y)
wire1 = ax.plot_wireframe(X,Y,u)
wire2 = ax.plot_wireframe(X,Y,v)
plt.show()
