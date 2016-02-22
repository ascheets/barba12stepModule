import numpy as np
import sympy
#from sympy import init_printing #include these when working with a notebook
#init_print(use_latex=True)
from sympy.utilities.lambdify import lambdify
import matplotlib.pyplot as plt


#setting up intial conditions, little more complicated this time...

#setting up symbolic variables
x, nu, t = sympy.symbols('x nu t')
phi = sympy.exp(-(x-4*t)**2/(4*nu*(t+1))) + sympy.exp(-(x-4*t-2*np.pi)**2/(4*nu*(t+1)))

phiprime = phi.diff(x)

u = -2*nu*(phiprime/phi)+4

ufunc = lambdify((t,x,nu),u) 

#domain declarations
nx = 101
nt = 100
dx = 2*np.pi/(nx-1)
nu = 0.07
dt = dx*nu

x = np.linspace(0, 2*np.pi, nx)
#u = np.empty(nx)
un = np.empty(nx)
t = 0

u = np.asarray([ufunc(t,x0,nu) for x0 in x])

plt.figure(figsize=(11,7), dpi=100)
plt.plot(x,u, marker='o', lw=2)
plt.xlim([0,2*np.pi])
plt.ylim([0,10])
plt.show()

for n in range(nt):
    un = u.copy()
    for i in range(nx-1):
        u[i] = un[i] -un[i] * dt/dx * (un[i] - un[i-1]) + nu*dt/dx**2*(un[i+1]-2*un[i]+un[i-1])
    u[-1] = un[-1] - un[-1] * dt/dx * (un[-1] - un[-2]) + nu*dt/dx**2*(un[0]-2*un[-1]+un[-2])

u_analytical = np.asarray([ufunc(nt*dt, xi, nu) for xi in x])

plt.figure(figsize=(11,7), dpi=100)
plt.plot(x,u,marker='o',lw=2,label='Computational')
plt.plot(x,u_analytical,label='Analytical')
plt.xlim([0,2*np.pi])
plt.ylim([0,10])
plt.legend()
plt.show()
