import numpy             as np
import numpy.random      as npr
import scipy.special     as sps
from matplotlib.pyplot   import figure
import matplotlib.pyplot as plt
import matplotlib        as mpl

sapin = (0.0353, 0.3216, 0.1569)

def signal_display(time_t,nsignal):

    plt.plot(time_t,nsignal,color = sapin, linewidth = 1);
    plt.grid()
    plt.xlabel('$t$ (s)')
    plt.tight_layout()

def planar_display(Kz, zt, zp, Nt = 500, Np = 500):
    
    thetas = np.linspace(1e-10,np.pi,Nt)
    phis   = np.linspace(0,2*np.pi,Np)
    
    # Display the spectrogram and its zeros
    figure(figsize = (8,6));
    plt.pcolormesh(phis, thetas,np.abs(Kz), shading = 'gouraud',cmap='gray')
    plt.scatter(zp,zt,s = 15,color = 'mistyrose')
    plt.xlabel(r'$\varphi$',fontsize = 45)
    plt.ylabel(r'$\vartheta$',fontsize = 45)
    # plt.xticks(np.pi/3*np.arange(7), ["$0$", "$\pi/3$", "$2\pi/3$", "$\pi$", "$4\pi/3$", "$5\pi/3$", "$2\pi$"]);
    # plt.yticks(np.pi/6*np.arange(7), ["$0$", "$\pi/6$", "$\pi/3$", "$\pi/2$", "$2\pi/3$", "$5\pi/6$", "$\pi$"]);
    plt.xticks(np.pi*np.arange(3), ["$0$", "$\pi$", "$2\pi$"]);
    plt.yticks(np.pi/2*np.arange(3), ["$0$", "$\pi/2$", "$\pi$"]);
    plt.tight_layout()

    
def spherical_display(Kz, zt, zp, Nt = 500, Np = 500):

    # maps of angles 
    thetas = np.linspace(1e-10,np.pi,Nt)
    phis   = np.linspace(0,2*np.pi,Np)
    Thetas, Phis = np.meshgrid(thetas,phis,indexing='ij')
                      
    # Cartesian coordinates of the unit sphere
    XX         = np.sin(Thetas) * np.cos(Phis)
    YY         = np.sin(Thetas) * np.sin(Phis)
    ZZ         = np.cos(Thetas)

    # Kravchuk transform of white noise
    test       = np.abs(Kz)
    tmin, tmax = test.min(), test.max()
    test       = (test - tmin)/(tmax - tmin)

    # Plot the surface
    fig = plt.figure(figsize=4*plt.figaspect(1.))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(XX, YY, ZZ, rstride=1, cstride=1, zorder=-1, facecolors=mpl.cm.gray(test), shade = False);
    ax.set_axis_off();  # Turn off the axis planes

    elev = -60
    azim = -145
    ax.view_init(elev, azim)

    # Add zeros
    xs = []
    ys = []
    zs = []
    eps=1e-3 # used to shift zeros towards the viewer, and bypass a rendering bug on zorder

    for i in range(len(zt)):

        # spherical coordinates
        phi   = zp[i]
        theta = zt[i]

        # Cartesian coordinates
        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)

        for j in range(2,3):
            npr.seed(j)
            n = npr.randn(3)
            if np.dot(n, np.array([x,y,z]))>=0: # that is, if the zero is on the hemisphere we can see
                #            if x-y+z>=0: # that is, if the zero is on the hemisphere we can see
                t = np.array([x,y,z]) + eps*n
                xs.append(t[0]) # shift the zeros towards the viewer
                ys.append(t[1])
                zs.append(t[2])

        ax.scatter(xs, ys, zs, s=75, zorder=-1, color="mistyrose")
    plt.tight_layout()