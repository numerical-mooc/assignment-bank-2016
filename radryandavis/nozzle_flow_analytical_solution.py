import numpy
from scipy import optimize

def f(M, A):
    '''Area/Mach relation in a form suitable for the bisection method solver.'''
    return A**2 - (1/M**2) * ((2/(1.4+1)) * (1 + ((1.4-1)/2) * M**2))**((1.4+1)/(1.4-1))

def compute_M(x, A):
    '''Computes Mach number in the nozzle.
    
    Parameters:
    ----------
    x : 1D array of float
        x grid
    A : 1D array of float
        nozzle area at each grid point
    
    Returns:
    -------
    M : 1D array of float
        Mach number at each grid point
    '''
    
    M = numpy.zeros_like(x)
    
    for i, x_val in enumerate(x):
        if x_val < 1.5:
            M[i] = optimize.bisect(f, 1e-6, 1, args=(A[i],))
        elif x_val == 1.5:
            M[i] = 1
        else:
            M[i] = optimize.bisect(f, 1, 50, args=(A[i],))
    
    return M

def nozzle_flow_analytical_solution(x, A):
    '''Computes the analytical solution to nozzle flow.
    
    Parameters:
    ----------
    x : 1D array of float
        x grid
    A : 1D array of float
        nozzle area at each grid point
    
    Returns:
    -------
    M : 1D array of float
        Mach number at each grid point
    p : 1D array of float
        pressure at each grid point
    rho : 1D array of float
          density at each grid point
    T : 1D array of float
        temperature at each grid point
    '''
    
    M = compute_M(x, A)
    p = (1 + ((1.4-1)/2) * M**2)**(-1.4/(1.4-1))
    rho = (1 + ((1.4-1)/2) * M**2)**(-1/(1.4-1))
    T = (1 + ((1.4-1)/2) * M**2)**(-1)
    
    return M, p, rho, T
