#step 12
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

nx = 41
ny = 41
nt = 10
nit=50
c = 1
dx = 2.0/(nx-1)
dy = 2.0/(ny-1)
x = np.linspace(0,2,nx)
y = np.linspace(0,2,ny)
X,Y = np.meshgrid(x,y)

rho = 1
nu = .1
F = 1
dt = .001

u = np.zeros((ny, nx))
#u[ny/4:ny/2,nx/4:nx/2] = .01
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

    #periodic BC pressure @ x = 2,0 so...
    #explicitly calculate b @ x = 2
    b[1:-1,-1]=(rho*(1/dt*\
                      ((u[1:-1,0]-u[1:-1,-2])/(2*dx)+ #gradu,gradx
                       (v[2:,-1]-v[0:-2,-1])/(2*dy))- #gradv,grady
                      ((u[1:-1,0]-u[1:-1,-2])/(2*dx))**2- #gradu,gradx squared
                      2*((u[2:,-1]-u[0:-2,-1])/(2*dy)*(v[1:-1,0]-v[1:-1,-2])/(2*dx))- #2*(gradu,grady)(gradv,gradx)
                      ((v[2:,-1]-v[0:-2,-1])/(2*dy))**2)) #gradv,grady squared

    #explicitly calculate b @ x = 0
    b[1:-1,0]=(rho*(1/dt*\
                      ((u[1:-1,1]-u[1:-1,-1])/(2*dx)+ #gradu,gradx
                       (v[2:,0]-v[0:-2,0])/(2*dy))- #gradv,grady
                      ((u[1:-1,1]-u[1:-1,-1])/(2*dx))**2- #gradu,gradx squared
                      2*((u[2:,0]-u[0:-2,0])/(2*dy)*(v[1:-1,1]-v[1:-1,-1])/(2*dx))- #2*(gradu,grady)(gradv,gradx)
                      ((v[2:,0]-v[0:-2,0])/(2*dy))**2)) #gradv,grady squared



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
        
        #explicity calculate pressure at x = 2,0 (periodic BCS)
        p[1:-1,-1] =( ((pn[1:-1,0]+pn[1:-1,-2])*dy**2+(pn[2:,-1]+pn[0:-2,-1])*dx**2)/(2*(dx**2+dy**2)) - #remnants of laplacian pressure business
                        rho*dx**2*dy**2/(2*(dx**2+dy**2))* #weird term that comes about from transposing
                        b[1:-1,-1] ) #source term

        p[1:-1,0] =( ((pn[1:-1,1]+pn[1:-1,-1])*dy**2+(pn[2:,0]+pn[0:-2,0])*dx**2)/(2*(dx**2+dy**2)) - #remnants of laplacian pressure business
                        rho*dx**2*dy**2/(2*(dx**2+dy**2))* #weird term that comes about from transposing
                        b[1:-1,0] ) #source term

        #wall boundary conditions, pressure
        p[-1,:] =p[-2,:] ##dp/dy = 0 at y = 2
        p[0,:] = p[1,:]  ##dp/dy = 0 at y = 0
        
        #previous conditions along x boundaries dealt with above
        #p[:,0]=p[:,1]    ##dp/dx = 0 at x = 0
        #p[:,-1]=0        ##p = 0 at x = 2
        
    return p

def channelFlow(nt, u, v, dt, dx, dy, p, rho, nu):
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
                            dt/dy**2*(un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1]))+ #two diffusive terms
                        F*dt) #source term

        v[1:-1,1:-1] =( vn[1:-1,1:-1]- #from FD in time
                        un[1:-1,1:-1]*dt/dx*(vn[1:-1,1:-1]-vn[1:-1,0:-2])- #gradv,gradx
                        vn[1:-1,1:-1]*dt/dy*(vn[1:-1,1:-1]-vn[0:-2,1:-1])- #gradv,grady
                        dt/(2*rho*dy)*(p[2:,1:-1]-p[0:-2,1:-1])+ #pressure term, gradp,grady
                        nu*(dt/dx**2*(vn[1:-1,2:]-2*vn[1:-1,1:-1]+vn[1:-1,0:-2])+ 
                            (dt/dy**2*(vn[2:,1:-1]-2*vn[1:-1,1:-1]+vn[0:-2,1:-1]))) ) #two diffusive terms

        #holy boundary conditions...

        #periodic BCs for u,v @ x = 2,0

        ####Periodic BC u @ x = 2     
        u[1:-1,-1] = un[1:-1,-1]-\
                     un[1:-1,-1]*dt/dx*(un[1:-1,-1]-un[1:-1,-2])-\
                     vn[1:-1,-1]*dt/dy*(un[1:-1,-1]-un[0:-2,-1])-\
                     dt/(2*rho*dx)*(p[1:-1,0]-p[1:-1,-2])+\
                     nu*(dt/dx**2*(un[1:-1,0]-2*un[1:-1,-1]+un[1:-1,-2])+\
                         dt/dy**2*(un[2:,-1]-2*un[1:-1,-1]+un[0:-2,-1]))+F*dt

        ####Periodic BC u @ x = 0
        u[1:-1,0] = un[1:-1,0]-\
                    un[1:-1,0]*dt/dx*(un[1:-1,0]-un[1:-1,-1])-\
                    vn[1:-1,0]*dt/dy*(un[1:-1,0]-un[0:-2,0])-\
                    dt/(2*rho*dx)*(p[1:-1,1]-p[1:-1,-1])+\
                    nu*(dt/dx**2*(un[1:-1,1]-2*un[1:-1,0]+un[1:-1,-1])+\
                        dt/dy**2*(un[2:,0]-2*un[1:-1,0]+un[0:-2,0]))+F*dt

        ####Periodic BC v @ x = 2
        v[1:-1,-1] = vn[1:-1,-1]-\
                     un[1:-1,-1]*dt/dx*(vn[1:-1,-1]-vn[1:-1,-2])-\
                     vn[1:-1,-1]*dt/dy*(vn[1:-1,-1]-vn[0:-2,-1])-\
                     dt/(2*rho*dy)*(p[2:,-1]-p[0:-2,-1])+\
                     nu*(dt/dx**2*(vn[1:-1,0]-2*vn[1:-1,-1]+vn[1:-1,-2])+\
                         (dt/dy**2*(vn[2:,-1]-2*vn[1:-1,-1]+vn[0:-2,-1])))

        ####Periodic BC v @ x = 0
        v[1:-1,0] = vn[1:-1,0]-\
                    un[1:-1,0]*dt/dx*(vn[1:-1,0]-vn[1:-1,-1])-\
                    vn[1:-1,0]*dt/dy*(vn[1:-1,0]-vn[0:-2,0])-\
                    dt/(2*rho*dy)*(p[2:,0]-p[0:-2,0])+\
                    nu*(dt/dx**2*(vn[1:-1,1]-2*vn[1:-1,0]+vn[1:-1,-1])+\
                        (dt/dy**2*(vn[2:,0]-2*vn[1:-1,0]+vn[0:-2,0])))   

        #"Wall BC" u,v = 0 at y = 0,2 -->no slip boundary condition
        u[0,:] = 0
        #u[:,0] = 0 
        #u[:,-1] = 0
        u[-1,:] = 0    #set velocity on upper, lower walls = 0
        v[0,:] = 0
        v[-1,:]=0
        #v[:,0] = 0
        #v[:,-1] = 0 #x boundary conditions periodic, dealt with above...
        
        
    return u, v, p

def channel(nt):

    u = np.zeros((ny, nx))
    v = np.zeros((ny, nx))
    p = np.zeros((ny, nx))
    b = np.zeros((ny, nx))
    u, v, p = channelFlow(nt, u, v, dt, dx, dy, p, rho, nu)
    fig = plt.figure(figsize=(11,7), dpi=100)
    #plt.contourf(X,Y,p,alpha=0.5)    ###plotting the pressure field as a contour
    #plt.colorbar()
    #plt.contour(X,Y,p)               ###plotting the pressure field outlines
    plt.quiver(X[::3,::3],Y[::3,::3],u[::3,::3],v[::3,::3]) ##plotting velocity ::3 --> only plot every 3rd data point
    plt.xlabel('X')
    plt.ylabel('Y')


    return

channel(25)
channel(50)
channel(75)
channel(100)
channel(499)
channel(1500)
plt.show()
