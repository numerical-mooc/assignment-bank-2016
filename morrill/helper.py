"""
helper.py
Contains functions to be used in Pressure_waves_in_a_box.ipynb
the final project in MAE 6286

Author: Michael Orrill

Last edit: 12/13/16
######################################################################
# Contains the following functions:                                  #
#                                                                    #
# get_coeff() --- calculates coefficients for Fourier series         #
# get_Fseries() - calculates Fourier series                          #
# animate1d() --- function for FuncAnimate for contour animation     #
# animate13d() -- function for FuncAnimate for 1 surface plot        #
# animate33d() -- function for FuncAnimate for 3 surface plots       #
# animate53d() -- function for FuncAnimate for 5 surface plots       #
# plot1d() ------ generates figure and image object for contour plot #
# plot13d() ----- generates figure and axis for 1 surface plot       # 
# plot33d() ----- generates figure and axis for 3 surface plots      #
# plot53d() ----- generates figure and axis for 5 surface plot2      # 
######################################################################
"""
import numpy
from matplotlib import pyplot, cm, rcParams
from mpl_toolkits.mplot3d import Axes3D
rcParams['font.size'] = 14
from matplotlib import animation
from IPython.display import HTML
# functions for fourier series:
def get_coeff(n,m,Lx,Ly,simple1=False,simple2=False):
    """Generates coefficient for Fourier series solution of 2D wave equation
    
    Parameters:
    -----------
    n:  int
         x term number in Fourier series
    m:  int
         y term number in Fourier series
    Lx: float
         Length in x direction
    Ly: float
         Length in y direction
        
    Returns:
    --------
    b: float
        coefficient to the Fourier series for term n,m
    """
    # for first simple comparison, f(x,y) = cos(x)cos(y)
    # http://www.wolframalpha.com/input/?i=int_0%5EL+cos(x)cos(b%2FL*x)dx
    if simple1==True:
        # things to simplify expresion
        pi = numpy.pi
        ax = n*pi
        ay = m*pi

        # for f(x,y) = cos(x)cos(y)
        bx = Lx*(Lx*numpy.cos(ax)*numpy.sin(Lx) - ax*numpy.sin(ax)*numpy.cos(Lx))/(Lx**2 - ax**2)
        by = Ly*(Ly*numpy.cos(ay)*numpy.sin(Ly) - ay*numpy.sin(ay)*numpy.cos(Ly))/(Ly**2 - ay**2)
        b = (4/(Lx*Ly))*bx*by

    # for second simple comparison f(x,y) = 2*sin(pi/4 x)sin(pi/4 y)
    # http://www.wolframalpha.com/input/?i=int_0%5EL+2*sin(b%2F4*x)cos(b%2FL*x)dx
    elif simple2==True:
        # things to simplify expresion
        pi = numpy.pi
        if m==1:
            m=0.99999
        if n==1:
            n=0.99999
        # for f(x,y) = 2sin(pi/4 x)sin(pi/4 y)
        bx = -4*Lx*(Lx*numpy.cos(pi*Lx/4)*numpy.sin(pi*n) - 4*n*numpy.sin(pi*Lx/4)*numpy.cos(pi*n))/(pi*(Lx**2 - 16*n**2))
        by = -4*Ly*(Ly*numpy.cos(pi*Ly/4)*numpy.sin(pi*m) - 4*m*numpy.sin(pi*Ly/4)*numpy.cos(pi*m))/(pi*(Ly**2 - 16*m**2))
        b = 2*(4/(Lx*Ly))*bx*by

    return b

def get_Fseries(x,y,n,m,Lx,Ly,c,nt,nx,ny,dt,f):
    """Generates Fourier series solution of 2D wave equation
    
    Parameters:
    -----------
    x:  2D array of float
         x mesh
    y:  2D array of float
         y mesh
    n:  int
         x term number in Fourier series
    m:  int
         y term number in Fourier series
    Lx: float
         Length in x direction
    Ly: float
         Length in y direction
    c:  float
         wavespeed
    nt: int
         number of time steps
    nx: int
         number of x direction steps
    ny: int
         number of y direction steps
    dt: float
         temporal discretization
    f:  int
         indecates which initial condition to use
        
    Returns:
    --------
    p: float
        Fourier series solution for 2D wave equation
    """
    # initiliaze p array
    p = numpy.ndarray((nt,ny,nx))
    if f==1:
        # first simple f(x,y) f(x,y) = cos(x)cos(y) and Neumann BCs
        for t in range(0,nt):
            for mi in range(0,m):
                for ni in range(0,n):
                    # things to simplify expressions:
                    ax = ni*numpy.pi/Lx
                    ay = mi*numpy.pi/Ly
                    at = numpy.sqrt(ax**2 + ay**2)
                
                    b = get_coeff(ni,mi,Lx,Ly,simple1=True) # coefficient of Fourier series
                    p[t] += b*numpy.cos(at*c*dt*t)*numpy.cos(ax*x)*numpy.cos(ay*y) # Fourier series summation
    elif f==2:
        # second simple f(x,y) = 2sin(pi/4 x)sin(pi/4 y) and Dirichlet BCs
        for t in range(0,nt):
            for mi in range(1,m+1):
                for ni in range(1,n+1):
                    # things to simplify expressions:
                    ax = ni*numpy.pi/Lx
                    ay = mi*numpy.pi/Ly
                    at = numpy.sqrt(ax**2 + ay**2)
                
                    b = get_coeff(ni,mi,Lx,Ly,simple2=True) # coefficient of Fourier series
                    p[t] += b*numpy.cos(at*c**2*dt*t)*numpy.sin(ax*x)*numpy.sin(ay*y) # Fourier series summation
    return p
	
# functions for animations
def animate1d(data):
    """Animation function for 1d contour plot
    
    Parameters:
    -----------
    data: numpy.ndarray of float
          contains solution to leapfrog at each time step
    Returns:
    --------
    im: pyplot.imshow object
    """
    im.set_array(data)
    return im,

def animate13d(i,p,lines):
    """Animation function for 1 3d surface plot
    
    Parameters:
    -----------
    i: int
        iterator
    p:    numpy.ndarray of float
          contains solution to leapfrog at each time step
    lines: the line objects for plotting
    Returns:
    --------
    line: object for plotting surface
    """
    ax.clear()
    line = ax.plot_surface(X,Y,p[i],cmap=cm.viridis)
    ax.set_zlim(0, 1)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    return line,

def animate33d(i,ps,pa,dif,lines):
    """animation function to comparing 3 3d surface plots
    
    Parameters:
    -----------
    i: int
        iterator
    ps: numpy.ndarray of float
           solution to the leapfrog method
    pa: numpy.ndarray of float
           solution to the analytical method
    dif: numpy.2D array of float
           difference between leap and analytical
           
    Returns:
    --------
    line: line object for plotting the surfaces
    """
    zmax = pa.max()
    zmin = pa.min()
    
    axs.clear()
    liness = axs.plot_surface(X,Y,ps[i],rstride=5,cstride=5,cmap=cm.viridis)
    axs.set_zlim(zmin,zmax)
    axs.set_xlabel('$x$')
    axs.set_ylabel('$y$')
    axs.set_title('Leapfrog')
    
    axa.clear()
    linea = axa.plot_surface(X,Y,pa[i],rstride=5,cstride=5,cmap=cm.viridis)
    axa.set_zlim(zmin,zmax)
    axa.set_xlabel('$x$')
    axa.set_ylabel('$y$')
    axa.set_title('Fourier series')
    
    axd.clear()
    lined = axd.plot_surface(X,Y,dif[i],rstride=5,cstride=5,cmap=cm.viridis)
    axd.set_zlim(zmin,zmax)
    axd.set_xlabel('$x$')
    axd.set_ylabel('$y$')
    axd.set_title('Difference of Leapfrog/Fourier')
    
    line = (liness,linea,lined)
    return line,

def animate53d(i,ps,pa,pf,difsa,diffa,lines):
    """animation function to comparing 3 3d surface plots
    
    Parameters:
    -----------
    i: int
        iterator
    ps: numpy.ndarray of float
           solution to the leapfrog method
    pa: numpy.ndarray of float
           solution to the analytical method
    pf: numpy.ndarray of float
           solution to the Fourier series
    dif: numpy.2D array of float
           difference between leap and analytical
           
    Returns:
    --------
    line: line object for plotting the surfaces
    """
    zmax = pa.max()
    zmin = pa.min()
    
    axs.clear()
    liness = axs.plot_surface(X,Y,ps[i],rstride=5,cstride=5,cmap=cm.viridis)
    axs.set_zlim(zmin, zmax)
    axs.set_xlabel('$x$')
    axs.set_ylabel('$y$')
    axs.set_title('Leapfrog')
    
    axa.clear()
    linea = axa.plot_surface(X,Y,pa[i],rstride=5,cstride=5,cmap=cm.viridis)
    axa.set_zlim(zmin,zmax)
    axa.set_xlabel('$x$')
    axa.set_ylabel('$y$')
    axa.set_title('Analytical')
    
    axf.clear()
    linef = axf.plot_surface(X,Y,pf[i],rstride=5,cstride=5,cmap=cm.viridis)
    axf.set_zlim(zmin,zmax)
    axf.set_xlabel('$x$')
    axf.set_ylabel('$y$')
    axf.set_title('Fourier series')
    
    axda.clear()
    lineda = axda.plot_surface(X,Y,difsa[i],rstride=5,cstride=5,cmap=cm.viridis)
    axda.set_zlim(zmin,zmax)
    axda.set_xlabel('$x$')
    axda.set_ylabel('$y$')
    axda.set_title('difference of Leapfrog/Analytical')
    
    axdf.clear()
    linedf = axdf.plot_surface(X,Y,diffa[i],rstride=5,cstride=5,cmap=cm.viridis)
    axdf.set_zlim(zmin,zmax)
    axdf.set_xlabel('$x$')
    axdf.set_ylabel('$y$')
    axdf.set_title('difference of Fourier/Analytical')
    
    line = (liness,linea,linef,lineda,linedf)
    return line,

def plot1d(p):
    """Generates figure and image for plotting animations of contour plot
    
    Parameters:
    -----------
    p: numpy.ndarray of float

    Returns:
    --------
    fig: figure for plotting
    
    im: image object for plotting the contour
    
    """
    fig = pyplot.figure(figsize=(10,4))
    pyplot.xticks([]), pyplot.yticks([]);
    im = pyplot.imshow(p[0], cmap=cm.viridis) 
    return fig, im
    
def plot13d(p,X,Y):
    """Generates figure and image for plotting animations of 1 3d surface plot
    
    Parameters:
    -----------
    p: numpy.ndarray of float
        solution to leapfrog at each time step
    X: numpy.array of float
        spatial discretization in x
    Y: numpy.array of float
        spatial discretization in y
    Returns:
    --------
    fig: figure for plotting
    
    line: line object for plotting the surface
    
    ax: axis objects for plot
    
    writer: the artist to draw the plot
    
    """
    fig = pyplot.figure(figsize=(10,4))
    ax = fig.add_subplot(111,projection='3d')

    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=25, metadata=dict(artist='Me'), bitrate=2800)

    line = ax.plot_surface(X,Y,p[0],cmap=cm.viridis)
    ax.set_zlim(0, 1)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    
    return fig, line, ax, writer

def plot33d(leap,analytical,dif,X,Y):
    """Generates figure and axi for plotting animations of 3 3d plots
    
    Parameters:
    -----------
    leap: numpy.ndarray of float
           solution to the leapfrog method
    analytical: numpy.ndarray of float
           solution to the analytical method
    dif: numpy.2D array of float
           difference between leap and analytical
    X: numpy.array of float
        spatial discretization in x
    Y: numpy.array of float
        spatial discretization in y      
    Returns:
    --------
    fig: figure for plotting
    
    line: line object for plotting the surfaces
    
    axs, axa, axd: axis objects for all 3 plots
    
    writer: the artist to draw the plots
    """
    zmax = analytical.max()
    zmin = analytical.min()
    
    fig = pyplot.figure(figsize=(10,10))
    axs = fig.add_subplot(221,projection='3d')
    axa = fig.add_subplot(222,projection='3d')
    axd = fig.add_subplot(223,projection='3d')
    pyplot.tight_layout()

    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=25, metadata=dict(artist='Me'), bitrate=2800)

    line = (axs.plot_surface(X,Y,leap[0],rstride=5,cstride=5,cmap=cm.viridis),\
            axa.plot_surface(X,Y,analytical[0],rstride=5,cstride=5,cmap=cm.viridis),\
            axd.plot_surface(X,Y,dif[0],rstride=5,cstride=5,cmap=cm.viridis))
    axs.set_zlim(zmin, zmax)
    axa.set_zlim(zmin, zmax)
    axd.set_zlim(zmin, zmax)

    return fig, line, axs, axa, axd, writer
    
def plot53d(leap,analytical,fourier,difla,diffa,X,Y):
    """Generates figure and axi for plotting animations of 3 3d plots
    
    Parameters:
    -----------
    leap: numpy.ndarray of float
           solution to the leapfrog method
    analytical: numpy.ndarray of float
           solution to the analytical method
    fourier: numpy.ndarray of float
           fourier series solution
    difla: numpy.2D array of float
           difference between leap and analytical
    diffa: numpy.2D array of float
           difference between fourier and analytical      
    Returns:
    --------
    fig: figure for plotting
    
    line: line object for plotting the surfaces
    
    axs, axa, axd: axis objects for all 3 plots
    
    writer: the artist to draw the plots
    """
    xmin = X.min()
    xmax = X.max()
    ymin = Y.min()
    ymax = Y.max()
    zmin = analytical.min()
    zmax = analytical.max()
    
    fig = pyplot.figure(figsize=(16,10))
    axs = fig.add_subplot(231,projection='3d')
    axs.set_xlim(xmin,xmax)
    axs.set_ylim(ymin,ymax)
    axs.set_zlim(zmin,zmax)
    
    axa = fig.add_subplot(232,projection='3d')
    axa.set_xlim(xmin,xmax)
    axa.set_ylim(ymin,ymax)
    axa.set_zlim(zmin,zmax)
    
    axf = fig.add_subplot(233,projection='3d')
    axf.set_xlim(xmin,xmax)
    axf.set_ylim(ymin,ymax)
    axf.set_zlim(zmin,zmax)
    
    axda = fig.add_subplot(223,projection='3d')
    axda.set_xlim(xmin,xmax)
    axda.set_ylim(ymin,ymax)
    axda.set_zlim(zmin,zmax)
    
    axdf = fig.add_subplot(224,projection='3d')
    axdf.set_xlim(xmin,xmax)
    axdf.set_ylim(ymin,ymax)
    axdf.set_zlim(zmin,zmax)
    #pyplot.tight_layout()

    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=25, metadata=dict(artist='Me'), bitrate=2800)

    line = (axs.plot_surface(X,Y,leap[0],rstride=5,cstride=5,cmap=cm.viridis),\
            axa.plot_surface(X,Y,analytical[0],rstride=5,cstride=5,cmap=cm.viridis),\
            axf.plot_surface(X,Y,fourier[0],rstride=5,cstride=5,cmap=cm.viridis),\
            axda.plot_surface(X,Y,difla[0],rstride=5,cstride=5,cmap=cm.viridis),\
            axdf.plot_surface(X,Y,diffa[0],rstride=5,cstride=5,cmap=cm.viridis))

    return fig, line, axs, axa, axf, axda, axdf, writer