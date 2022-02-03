import numpy              as np
from scipy.io             import savemat, loadmat
import matplotlib.pyplot  as plt
from spherical_statistics import empirical_F, empirical_K

def the_test_statistic(S,S0,the_norm=2,r2 = int(1e4)):
    
    """
    test statistics with possible choice of the norm
    """
    
    # range of rs of interest
    ind_rs = np.arange(0,r2)
    
    # test statistics
    T = np.linalg.norm(S[ind_rs]-S0[ind_rs],ord=the_norm)
    
    return T
    
    
def the_test(signal,alpha,m,folder='samples',functional='F'):
    
    """
    summary statistics used to perform the test
    """
    
    #signal: observed data to be tested
    #folder: folder in which the white 
    #functional: either 'F' for the empty space function or 'K' for Ripley's K function
    
    k = int(alpha*(m+1))
    
    #Kravchuk analysis and functional statistics of the data
    if functional     == 'F':
        S,_           = empirical_F(signal)
    else:
        S,_           = empirical_K(signal)
    nr        = len(S)
    S_m       = np.zeros((m+1,nr))
    S_m[-1,:] = S
    
    #load the functional statistics of the noise samples
    for n in np.arange(m):
        ndict         = loadmat('../'+folder+'/noise_'+str(n)+'.mat',squeeze_me=True);
        if functional == 'F':
            S_m[n,:]  = ndict["Fn"]
        else:
            S_m[n,:]  = ndict["RipK"]
    
    # compute the summary statistics
    S0        = np.sum(S_m, axis = 0)/(m+1)
    T_m       = np.zeros(m)
    for i in np.arange(m):
            T_m[i]    = the_test_statistic(S_m[i,:],S0)
    T_m.sort()
    t_exp             = the_test_statistic(S,S0)
    
    if t_exp >= T_m[-k]:
        print('The null hypothesis is rejected. A signal is detected.')
    else:
        print('The null hypothesis cannot be rejected. No signal detected.')
    
    print(' ')
    print('--------------------------------------------------------------')
    print(' ')
    print('Threshold of rejection: %.2d' % (T_m[-k]))
    print('Value of the summary statistics: %.2d' % t_exp)
    return 
    