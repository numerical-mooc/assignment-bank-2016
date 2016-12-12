import numpy
def get_coeff(n,m,Lx,Ly):
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
       
    Note:
    Expression calculated here: 
    http://www.wolframalpha.com/input/?i=int_0%5EL+exp(-(x-0.5*L)%5E2%2F(2*L%5E2))cos(b%2FL*x)dx
    """

    # things to simplify expresion
    pi = numpy.pi
    ax = n*pi
    ay = m*pi

    # for f(x,y) = cos(x)cos(y)
    bx = Lx*(Lx*numpy.cos(ax)*numpy.sin(Lx) - ax*numpy.sin(ax)*numpy.cos(Lx))/(Lx**2 - ax**2)
    by = Ly*(Ly*numpy.cos(ay)*numpy.sin(Ly) - ay*numpy.sin(ay)*numpy.cos(Ly))/(Ly**2 - ay**2)
    b = (4/(Lx*Ly))*bx*by
    return b

def get_Fseries(x,y,n,m,Lx,Ly,c,nt,nx,ny,dt):
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
        
    Returns:
    --------
    p: float
        Fourier series solution for 2D wave equation
    """
    # initiliaze p array
    p = numpy.ndarray((nt,ny,nx))
    
    for t in range(0,nt):
        for mi in range(0,m):
            for ni in range(0,n):
                # things to simplify expressions:
                ax = ni*numpy.pi/Lx
                ay = mi*numpy.pi/Ly
                at = numpy.sqrt(ax**2 + ay**2)
                
                b = get_coeff(ni,mi,Lx,Ly) # coefficient of Fourier series
                p[t] += b*numpy.cos(at*c*dt*t)*numpy.cos(ax*x)*numpy.cos(ay*y) # Fourier series summation
    return p