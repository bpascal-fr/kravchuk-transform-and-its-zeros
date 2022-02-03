import numpy             as np
import matplotlib.pyplot as plt
from kravchuk_transform  import the_transform, the_zeros

sapin = (0.0353, 0.3216, 0.1569)

def the_distance(theta1,phi1,theta2,phi2):
    
    """
    distance between the points of spherical coordinates (theta1,phi1) and (theta2,phi2)
    """
    
    s_xy = np.arccos(np.sin(theta1)*np.sin(theta2)*np.cos(phi1-phi2)+np.cos(theta1)*np.cos(theta2))
    return s_xy


def the_F_statistics(zt,zp,rs):
    
    """
    unbiased estimator of the empty space function for a point process on the sphere
    """
    
    # number of zeros
    Nz = zt.shape[0]
    
    # the grid for the quadrature of the uniform measure on the sphere
    N_theta  = 4*int(np.sqrt(Nz))
    N_phi    = 4*int(np.sqrt(Nz))
    N_u      = N_theta*N_phi   
    gthetas  = np.arccos(np.linspace(1,-1,N_theta))
    gphis    = np.linspace(0,2*np.pi,N_phi)

    
    # distances from each reference point to each zero of the Kravchuk transform
    S_xu = np.zeros((Nz,N_u))
    for nz in np.arange(Nz):
        nu   = 0 
        for gtheta in gthetas:
            for gphi in gphis:
                S_xu[nz,nu] = the_distance(zt[nz],zp[nz],gtheta,gphi)
                nu  += 1

    # distance to the nearest point of the point process
    inf_xu = np.amin(S_xu,axis = 0)
    inf_xu = np.reshape(inf_xu,(N_u,))
    
    # number of points at distance less than r
    Fr     = np.zeros(np.shape(rs))
    ir     = 0
    for r in rs:
        Fr [ir] = np.sum(inf_xu < r)/(N_u)
        ir      += 1 
    
    return Fr, rs


def empirical_F(signal,disp = False):
    
    Ts       = the_transform(signal)
    N        = len(signal)-1
    zt, zp   = the_zeros(Ts,N)
    r        = np.linspace(0,2*np.pi/np.sqrt(N),10**4)
    F,_      = the_F_statistics(zt,zp,r)
    
    if disp:
        plt.plot(r,F,color = sapin);
        plt.grid()
        plt.xlabel('$r$')
        plt.xlabel('$F(r)$')
        plt.tight_layout()
        
    return F, r

def the_K_statistics(zt,zp,rs=np.linspace(0,np.pi,10**4)):
    
    """
    unbiased estimator of Ripley's K function of a point process on the sphere
    """
    
    # number of zeros
    Nz = zt.shape[0]
    
    
    # distances between pairs
    Npairs     = int(Nz*(Nz-1)/2)
    S_pairs    = np.zeros(Npairs)
    n_pairs    = 0
    for nz in np.arange(Nz):
        for mz in np.arange(0,nz):
            S_pairs[n_pairs] = the_distance(zt[nz],zp[nz],zt[mz],zp[mz])
            n_pairs          += 1
    
    Ks         = np.zeros(len(rs))
    ir         = 0
    for r in rs:
        Ks[ir] = 2*(4*np.pi)**2*np.sum(S_pairs < r)/Nz**2
        ir     +=1
        
    return Ks, rs

def empirical_K(signal,disp = False):
    
    Ts     = the_transform(signal)
    N      = len(signal)-1
    zt, zp = the_zeros(Ts,N)
    K, r   = the_K_statistics(zt,zp)
    
    if disp:
        plt.plot(r,K,color = sapin);
        plt.grid()
        plt.xlabel('$r$')
        plt.xlabel('$K(r)$')
        plt.tight_layout()
        
    return K, r