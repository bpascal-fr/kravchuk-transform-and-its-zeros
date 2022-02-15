import numpy               as np
import numpy.random        as npr
import os
from scipy.io                        import savemat, loadmat

from chirp_signals          import the_chirp, the_white_noise, the_noisy_chirp
from kravchuk_transform     import the_transform, the_zeros
from kravchuk_display       import signal_display, planar_display, spherical_display
from spherical_statistics   import the_distance, the_F_statistics, the_K_statistics, empirical_F, empirical_K

from stft_transform  import the_stft_transform, the_stft_zeros, stft_display
from spatialstats    import spatialStatsFromR

def noise_samples(N,alpha = 0.05,m=199, time_t=[],folder = 'samples'):
    
    """
    generate m samples of complex white Gaussian noise and perform all required analysis for zero-based detection
    """
    
    #N+1: number of points in the white noise
    #alpha: desired significance level of the test
    #m:   number of realizations for the Monte Carlo envelope test
    #folder: name of the folder in which the white noise samples are stored
    
    # create a folder to save the white noise samples
    if os.path.isdir('../'+folder) == False:
        os.mkdir('../'+folder)
    if len(time_t) == 0:
        time_t     = np.arange(N+1)
        
    for n in np.arange(m):

        # sample white Gaussian noise
        wnoise   = the_white_noise(N)
        
        # Kravchuk spectrogram analysis
        # compute its Kravchuk transform
        Kn       = the_transform(wnoise)
        # find the zeros
        znt, znp = the_zeros(Kn,N)
        # explored range of r
        rs       = np.linspace(0,2*np.pi/np.sqrt(N),10**4)
        # compute the F-function
        Fn,_     = the_F_statistics(znt,znp,rs)
        # compute Ripley's K function
        RipK,_   = the_K_statistics(znt,znp)
        
        # Standard time-frequency analysis
        # compute the stft transform
        Vn, fn   = the_stft_transform(wnoise, time_t)
        # find its zeros 
        vnt, vnf = the_stft_zeros(Vn,time_t,fn)
        # compute the spatial statistics of the STFT
        r_K, KR, r_F, FR = spatialStatsFromR(time_t,fn,vnt,vnf)

        # Save the data
        mdict = {"N" : N, "wnoise" : wnoise, "time_t" : time_t, "Kn" : Kn, "znt" : znt, "znp" : znp,
                "rs" : rs, "Fn" : Fn, "RipK" : RipK, "Vn" : Vn, "fn" : fn, "vnt" : vnt, "vnf" : vnf,
                "r_K" : r_K, "KR" : KR, "r_F" : r_F, "FR" : FR}
        savemat('../samples/noise_'+str(n)+'.mat',mdict)
        
    return alpha, m, folder
