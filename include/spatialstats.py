import numpy as np
import matplotlib.pyplot as plt

import rpy2.robjects as robjects
# Activate automatic conversion of numpy floats and arrays to corresponding R objects
from rpy2.robjects import numpy2ri
numpy2ri.activate() #numpy2ri.deactivate()

from spatstat_interface.utils import to_pandas_data_frame
from spatstat_interface.interface import SpatstatInterface

spatstat = SpatstatInterface(update=True)
spatstat.import_package("core", "geom", update=True)

def spatialStatsFromR(time_t,f,zt,zf):

    sample = np.array([zt, zf])

    tbound  = robjects.FloatVector([time_t[0], time_t[-1]])
    fbound  = robjects.FloatVector([f[0], f[-1]])
    window  = spatstat.geom.owin(xrange=tbound, yrange=fbound)

    X       = spatstat.geom.ppp(*sample, window=window)

    numpy2ri.deactivate()
    Kest_r  = spatstat.core.Kest(X)
    Kest_df = to_pandas_data_frame(Kest_r)


    r_K     = np.array(Kest_df["r"])
    KR      = np.array(Kest_df["border"])
    
    Fest_r  = spatstat.core.Fest(X)
    Fest_df = to_pandas_data_frame(Fest_r)


    r_F     = np.array(Fest_df["r"])
    FR      = np.array(Fest_df["km"])
    
    numpy2ri.activate()
    
    return r_K, KR, r_F, FR

