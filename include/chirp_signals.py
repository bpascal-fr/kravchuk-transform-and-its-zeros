import numpy             as np
import matplotlib.pyplot as plt
import numpy.random      as npr

sapin = (0.0353, 0.3216, 0.1569)

def the_chirp(N,observe = 20,duration = 15):
    
    """
    chirp of unit energy for the 2-norm
    """
    
    # N+1: number of samples
    # observe: total duration of the measurement
    # duration: duration of the phenomenon of interest
    
    
    time_t   = np.linspace(-observe, observe, N+1)
    envelop  = np.zeros(N+1)
    envelop[np.abs(time_t)<duration] = np.exp(-25/(duration**2-time_t[np.abs(time_t)<duration]**2))
    freq     = 0.5 + (time_t + duration)*(1 - 0.5)/observe
    chirp    = np.sin(2*np.pi*freq*time_t)
    signal   = chirp*envelop
    signal   /= np.linalg.norm(signal,ord=2)

        
    return signal, time_t

def the_white_noise(N):
    
    """
    complex white Gaussian noise of length N+1 normalized for the 2-norm
    """
    
    # N+1: number of samples
    
    wnoise = (np.random.randn(N+1)+1j*np.random.randn(N+1))/np.sqrt(2)
    wnoise /= np.linalg.norm(wnoise,ord=2)
    
    return wnoise

def the_noisy_chirp(N,snr = 1,observe=20,duration=15,disp = False):
    
    """
    noisy chirp
    """
    
    # N+1: number of samples
    # snr: signal-to-noise ratio
    # observe: total duration of the measurement
    # duration: duration of the phenomenon of interest
    
    signal, time_t = the_chirp(N,observe,duration)
    
    wnoise = the_white_noise(N)
    
    if snr > 0:
        nsignal = signal + 1/snr*wnoise
    else:
        nsignal = wnoise
        
    if disp:
        display_signal(nsignal,time_t)
        
    return nsignal, time_t

def display_signal(nsignal,time_t=np.array([])):
                   
    """
    display the real part of a signal along time
    """
    
    # nsignal: complex signal to be analyzed of length N+1
    # time_t: (optional) time range of observation (default: [0,1, ...,N])
    
    if len(time_t) == 0:
        time_t = np.arange(len(nsignal))
        
    plt.plot(time_t,nsignal.real,color = sapin, linewidth = 1);
    plt.grid()
    plt.xlabel('$t$ (s)')
    plt.tight_layout()
    
    return time_t