#optimizing loops with numba
from numba import autojit

@autojit
def laplace2d_numba(p, b, y, dx, dy, l1norm_target):
    l1norm = 1
    pn = np.empty_like(p)

    while l1norm > l1norm_target :
        pn = p.copy()

        #note, i --> cols, j ---> rows
        p[1:-1,1:-1] = (dx**2*(pn[2:,1:-1] + pn[0:-2,1:-1]) + \
                        dy**2*(pn[1:-1,2:] + pn[1:-1,0:-2]) - \
                        dx**2*dy**2*(b[1:-1,1:-1]))/(2*(dx**2 + dy**2))

        p[:,0] = 0
        p[:,-1] = 0
        p[0,:] = 0
        p[-1,:] = 0
        l1norm = (np.sum(np.abs(p[:])-np.abs(pn[:])))/np.sum(np.abs(pn[:]))

    return p
