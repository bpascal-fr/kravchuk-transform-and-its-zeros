import numpy             as np
import scipy.signal      as signal
import scipy.fft         as fft
import matplotlib.pyplot as plt

def the_stft_transform(x, time_t=[]):
    
    """
    standard Short-Time Fourier Transform with a circular Gaussian window
    """
    
    # x:      complex signal to be analyzed with N+1 samples
    # time_t: time_t: (optional) time range of observation (default: [0,1, ...,N])
    
    if len(time_t) == 0:
        time_t     = np.arange(len(x))
        
    dt = time_t[1]-time_t[0]

    # Gaussian window of unit energy
    Ng       = len(x)
    g        = signal.gaussian(Ng, np.sqrt((Ng)/2/np.pi))
    g        = g/g.sum()

    # Compute the Gaussian Short-Time Fourier Transform
    fx, _, Vx = signal.stft(x, fs = 1/dt, window=g, nperseg=Ng, noverlap=Ng-1, return_onesided=False)
    
    # remap the frequencies in a symmetric fashion
    fx = fft.fftshift(fx)
    Vx = fft.fftshift(Vx,axes = (0,))
    
    return Vx, fx


def the_stft_zeros(Vx,time_t,fx,feedback = False):

    # compute the zeros
    zx, zy    = extr2min(np.abs(Vx))
    
    # coordinates of the zeros on the sphre
    zt = time_t[zy]
    zf = fx[zx]
    
    # give some feedback
    if feedback:
        print('Local minima method has found '+str(zx.shape[0])+' zeros.')
    
    return zt, zf


def stft_display(Vx,time_t,fx,zt,zf):
    
    plt.pcolormesh(time_t,fx,np.abs(Vx), shading='gouraud', cmap = 'gray')
    plt.scatter(zt,zf,s = 15,color = 'mistyrose')
    plt.xlabel(r'$t$ (s)')
    plt.ylabel(r'$\omega$ (Hz)')
    plt.tight_layout()
    
    
def extr2min(M):
    
    """
    find spectrogram zeros by the Minimal Grid Neighbors method
    """
    
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
