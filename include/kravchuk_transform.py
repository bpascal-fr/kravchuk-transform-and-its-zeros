import numpy as np
import scipy.special as sps

def the_transform(x,Nt = 500,Np = 500):

    N  = x.shape[0] 
    thetas = np.linspace(1e-10,np.pi,Nt)
    phis   = np.linspace(0,2*np.pi,Np)
    
    # duplicate the signal
    X = np.kron(x,np.ones((Np,Nt,1))).T
    
    # maps of angles 
    Thetas, Phis = np.meshgrid(thetas,phis,indexing='ij')
    
    # map of phase space variable
    z = np.cos(Thetas/2)/np.sin(Thetas/2)*np.exp(1j*Phis)
    Z = np.kron(z,np.ones((N,1,1)))
    
    # indices
    ell = np.arange(N)
    L = np.kron(ell,np.ones((Np,Nt,1))).T
    
    # compute the transform using the generative function trick
    normZ  = np.sqrt(2*(1+np.abs(Z)**2))
    sumand = X*np.sqrt(sps.comb(N,L))*((1-Z)/(normZ))**L*((1+Z)/(normZ))**(N-L)
    Kz = np.sum(sumand,0)
    
    return Kz


def the_zeros(Kz, N, Nt = 500, Np = 500,feedback = False):

    thetas = np.linspace(1e-10,np.pi,Nt)
    phis   = np.linspace(0,2*np.pi,Np)
    
    # compute the zeros
    zx, zy    = extr2min(np.abs(Kz))
    
    # coordinates of the zeros on the sphre
    zp = phis[zy]
    zt = thetas[zx]
    
    # give some feedback
    if feedback:
        print('Local minima method has found '+str(zx.shape[0])+' zeros compared to the '+str(N)+' expected.')
    
    return zt, zp
    

def extr2min(M):
    central = M[1:-1, 1:-1]
    mask = np.full(central.shape, True, dtype=bool)
    sub_indices = (
        (np.s_[2:], np.s_[1:-1]),
        (np.s_[:-2], np.s_[1:-1]),
        (np.s_[1:-1], np.s_[2:]),
        (np.s_[1:-1], np.s_[:-2]),
        (np.s_[:-2], np.s_[:-2]),
        (np.s_[:-2], np.s_[2:]),
        (np.s_[2:], np.s_[2:]),
        (np.s_[2:], np.s_[:-2]),
    )
    for I, J in sub_indices:
        np.logical_and(mask, central <= M[I, J], out=mask, where=mask)

    x, y = np.where(mask)
    return x + 1, y + 1

    
    