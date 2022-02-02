import numpy             as np
from scipy.io            import savemat, loadmat
import matplotlib.pyplot as plt

def the_test_statistic(S,S0,the_norm=2,r2 = int(1e4)):
    
    """
    test statistics in either inftinity of l2 norms
    """
    
    # range of rs of interest
    ind_rs = np.arange(0,r2)
    
    # test statistics
    T = np.linalg.norm(S[ind_rs]-S0[ind_rs],ord=the_norm)
    
    return T
    
    
