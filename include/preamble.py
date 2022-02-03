import sys
sys.path.append("../include")

import matplotlib                    as mpl
import matplotlib.pyplot             as plt
import numpy.random                  as npr
import statsmodels.stats.proportion  as stms
from mpl_toolkits.axes_grid1         import make_axes_locatable
from scipy.io                        import savemat, loadmat

from chirp_signals          import the_chirp, the_white_noise, the_noisy_chirp, display_signal
from kravchuk_transform     import the_transform, the_zeros
from kravchuk_display       import signal_display, planar_display, spherical_display
from spherical_statistics   import the_distance, the_F_statistics, the_K_statistics, empirical_F, empirical_K
from detection_test         import the_test_statistic, the_test
from stft_transform  import the_stft_transform, the_stft_zeros, stft_display
from white_noises    import noise_samples


mpl.rcParams['xtick.labelsize'] = 30;
mpl.rcParams['ytick.labelsize'] = 30;
mpl.rcParams['axes.titlesize'] = 30;
plt.rc('axes', labelsize=35);
plt.rc('legend', fontsize=30);
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
mpl.rcParams['font.family'] = 'roman'

sapin  = (0.0353, 0.3216, 0.1569)
carmin = (0.7294,0.0392,0.0392)
bleu   = (0.2, 0.2, 0.7020)
rose   = (0.56, 0.004, 0.32)
marron = (0.51, 0.3, 0.11)
jaune  = (0.93, 0.69, 0.13)