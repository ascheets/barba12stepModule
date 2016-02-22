#step 11
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

nx = 41
ny = 41
nt = 500
nit=50
c = 1
dx = 2.0/(nx-1)
dy = 2.0/(ny-1)
x = np.linspace(0,2,nx)
y = np.linspace(0,2,ny)
X,Y = np.meshgrid(x,y)

rho = 1
nu = .1
dt = .001

u = np.zeros((ny, nx))
v = np.zeros((ny, nx))
p = np.zeros((ny, nx)) 
b = np.zeros((ny, nx))

#building up the source term for poisson pressure equation
def buildUpB(b, rho, dt, u, v, dx, dy):
    
    b[1:-1,1:-1]=(rho*(1/dt*\
                      ((u[1:-1,2:]-u[1:-1,0:-2])/(2*dx)+ #gradu,gradx
                       (v[2:,1:-1]-v[0:-2,1:-1])/(2*dy))- #gradv,grady
                      ((u[1:-1,2:]-u[1:-1,0:-2])/(2*dx))**2- #gradu,gradx squared
                      2*((u[2:,1:-1]-u[0:-2,1:-1])/(2*dy)*(v[1:-1,2:]-v[1:-1,0:-2])/(2*dx))- #2*(gradu,grady)(gradv,gradx)
                      ((v[2:,1:-1]-v[0:-2,1:-1])/(2*dy))**2)) #gradv,grady squared

    return b

#function which ensures continuity by using pressure constraint to enforce zero velocity divergence
#iterative function, relaxes pressure in artificial time
def presPoisson(p, dx, dy, b):
    pn = np.empty_like(p)
    pn = p.copy()
    
    for q in range(nit):
        pn = p.copy()
        p[1:-1,1:-1] =( ((pn[1:-1,2:]+pn[1:-1,0:-2])*dy**2+(pn[2:,1:-1]+pn[0:-2,1:-1])*dx**2)/(2*(dx**2+dy**2)) - #remnants of laplacian pressure business
                        rho*dx**2*dy**2/(2*(dx**2+dy**2))* #weird term that comes about from transposing
                        b[1:-1,1:-1] ) #source term

        p[-1,:] =p[-2,:] ##dp/dy = 0 at y = 2
        p[0,:] = p[1,:]  ##dp/dy = 0 at y = 0
        p[:,0]=p[:,1]    ##dp/dx = 0 at x = 0
        p[:,-1]=0        ##p = 0 at x = 2
        
    return p

def cavityFlow(nt, u, v, dt, dx, dy, p, rho, nu):
    un = np.empty_like(u)
    vn = np.empty_like(v)
    b = np.zeros((ny, nx))
    
    for n in range(nt):
        un = u.copy()
        vn = v.copy()
        
        b = buildUpB(b, rho, dt, u, v, dx, dy)
        p = presPoisson(p, dx, dy, b)
        
        u[1:-1,1:-1] =( un[1:-1,1:-1]- #from FD in time
                        un[1:-1,1:-1]*dt/dx*(un[1:-1,1:-1]-un[1:-1,0:-2])- #gradu,gradx
                        vn[1:-1,1:-1]*dt/dy*(un[1:-1,1:-1]-un[0:-2,1:-1])- #gradv,grady
                        dt/(2*rho*dx)*(p[1:-1,2:]-p[1:-1,0:-2])+ #pressure term, gradp,gradx
                        nu*(dt/dx**2*(un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2])+
                        dt/dy**2*(un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1])) ) #two diffusive terms

        v[1:-1,1:-1] =( vn[1:-1,1:-1]- #from FD in time
                        un[1:-1,1:-1]*dt/dx*(vn[1:-1,1:-1]-vn[1:-1,0:-2])- #gradv,gradx
                        vn[1:-1,1:-1]*dt/dy*(vn[1:-1,1:-1]-vn[0:-2,1:-1])- #gradv,grady
                        dt/(2*rho*dy)*(p[2:,1:-1]-p[0:-2,1:-1])+ #pressure term, gradp,grady
                        nu*(dt/dx**2*(vn[1:-1,2:]-2*vn[1:-1,1:-1]+vn[1:-1,0:-2])+ 
                        (dt/dy**2*(vn[2:,1:-1]-2*vn[1:-1,1:-1]+vn[0:-2,1:-1]))) ) #two diffusive terms

        u[0,:] = 0
        u[:,0] = 0
        u[:,-1] = 0
        u[-1,:] = 1    #set velocity on cavity lid equal to 1
        v[0,:] = 0
        v[-1,:]=0
        v[:,0] = 0
        v[:,-1] = 0
        
        
    return u, v, p

def cavity(nt):

    u = np.zeros((ny, nx))
    v = np.zeros((ny, nx))
    p = np.zeros((ny, nx))
    b = np.zeros((ny, nx))
    u, v, p = cavityFlow(nt, u, v, dt, dx, dy, p, rho, nu)
    fig = plt.figure(figsize=(11,7), dpi=100)
    plt.contourf(X,Y,p,alpha=0.5)    ###plotting the pressure field as a contour
    plt.colorbar()
    plt.contour(X,Y,p)               ###plotting the pressure field outlines
    plt.quiver(X[::2,::2],Y[::2,::2],u[::2,::2],v[::2,::2]) ##plotting velocity
    plt.xlabel('X')
    plt.ylabel('Y')


    return

cavity(25)
cavity(50)
cavity(75)
cavity(100)
#cavity(700)
#cavity(1500)
plt.show()
