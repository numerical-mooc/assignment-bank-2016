import numpy
from matplotlib import pyplot, rcParams
from scipy import optimize
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

def compute_initial_conditions(x):
    '''Computes and returns the initial conditions for density, velocity, and temperature.
       Also computes and returns the area of the nozzle.
    
    Parameters:
    ----------
    x : 1D array of float
        x grid
    
    Returns:
    -------
    A : 1D array of float
        nozzle area
    rho : 1D array of float
          initial density
    V : 1D array of float
        initial velocity
    T : 1D array of float
        initial temperature
    '''
    
    A = 1 + 2.2 * (x - 1.5)**2
    
    rho = 1 - 0.3146 * x
    
    T = 1 - 0.2314 * x
    
    V = (0.1 + 1.09 * x) * numpy.sqrt(T)
    
    return A, rho, V, T

def solve_nozzle_flow(nt, nx, dx, A, rho, V, T):
    '''Solves the nozzle flow problem.
    
    Parameters:
    ----------
    nt : int
         number of time steps
    nx : int
         number of grid points
    dx : float
         grid spacing
    A : 1D array of float
        area at each grid point
    rho : 1D array of float
          density at each grid point
    V : 1D array of float
        velocity at each grid point
    T : 1D array of float
        temperature at each grid point
    
    Returns:
    -------
    t_n : 1D array of float
          array of time-step size
    rho_n : 2D array of float
            density array at each time step
    V_n : 2D array of float
          velocity array at each time step
    T_n : 2D array of float
          temperature array at each time step
    '''
    
    C = .5  # Courant number
    
    rho_dot_predictor = numpy.zeros(nx)
    V_dot_predictor = numpy.zeros(nx)
    T_dot_predictor = numpy.zeros(nx)
    
    rho_dot_corrector = numpy.zeros(nx)
    V_dot_corrector = numpy.zeros(nx)
    T_dot_corrector = numpy.zeros(nx)
    
    rho_next = numpy.zeros(nx)
    V_next = numpy.zeros(nx)
    T_next = numpy.zeros(nx)
    
    rho_n = numpy.zeros((nt,nx))
    V_n = numpy.zeros((nt,nx))
    T_n = numpy.zeros((nt,nx))
    
    t_n = numpy.zeros(nt)
    
    for t in range(nt):
        
        dt = numpy.amin(C * (dx/(numpy.sqrt(T) + V)))
        
        # MacCormack's Method: Predictor Step
        rho_dot_predictor[:-1] = -rho[:-1] * ((V[1:] - V[:-1])/dx) -\
                                  rho[:-1] * V[:-1] * ((numpy.log(A[1:]) - numpy.log(A[:-1]))/dx) -\
                                  V[:-1] * ((rho[1:] - rho[:-1])/dx)
        
        V_dot_predictor[:-1] = -V[:-1] * ((V[1:] - V[:-1])/dx) -\
                                (1/1.4) * (((T[1:] - T[:-1])/dx) +\
                                           (T[:-1]/rho[:-1]) * ((rho[1:] - rho[:-1])/dx))
        
        T_dot_predictor[:-1] = -V[:-1] * ((T[1:] - T[:-1])/dx) -\
                                (1.4-1) * T[:-1] * (((V[1:] - V[:-1])/dx) +\
                                                    V[:-1] * ((numpy.log(A[1:]) - numpy.log(A[:-1]))/dx))
        
        rho_predicted = rho.copy()
        V_predicted = V.copy()
        T_predicted = T.copy()
        
        rho_predicted[:-1] = rho[:-1] + rho_dot_predictor[:-1] * dt
        
        V_predicted[:-1] = V[:-1] + V_dot_predictor[:-1] * dt
        
        T_predicted[:-1] = T[:-1] + T_dot_predictor[:-1] * dt
        
        
        # MacCormack's Method: Corrector Step
        rho_dot_corrector[1:] = -rho_predicted[1:] * ((V_predicted[1:] - V_predicted[:-1])/dx) -\
                                 rho_predicted[1:] * V_predicted[1:] *\
                                 ((numpy.log(A[1:]) - numpy.log(A[:-1]))/dx) -\
                                 V_predicted[1:] * ((rho_predicted[1:] - rho_predicted[:-1])/dx)
        
        V_dot_corrector[1:] = -V_predicted[1:] * ((V_predicted[1:] - V_predicted[:-1])/dx) -\
                               (1/1.4) * (((T_predicted[1:] - T_predicted[:-1])/dx) +\
                                          (T_predicted[1:]/rho_predicted[1:]) *\
                                          ((rho_predicted[1:] - rho_predicted[:-1])/dx))
        
        T_dot_corrector[1:] = -V_predicted[1:] * ((T_predicted[1:] - T_predicted[:-1])/dx) -\
                               (1.4-1) * T_predicted[1:] *\
                               (((V_predicted[1:] - V_predicted[:-1])/dx) +\
                                V_predicted[1:] * ((numpy.log(A[1:]) - numpy.log(A[:-1]))/dx))
        
        rho_dot_average = .5 * (rho_dot_predictor + rho_dot_corrector)
        
        V_dot_average = .5 * (V_dot_predictor + V_dot_corrector)
        
        T_dot_average = .5 * (T_dot_predictor + T_dot_corrector)
        
        rho_next[1:-1] = rho[1:-1] + rho_dot_average[1:-1] * dt
        
        V_next[1:-1] = V[1:-1] + V_dot_average[1:-1] * dt
        
        T_next[1:-1] = T[1:-1] + T_dot_average[1:-1] * dt
        
        
        # Boundary Conditions
        # inflow
        rho_next[0] = 1
        V_next[0] = 2 * V_next[1] - V_next[2]
        T_next[0] = 1
        
        # outflow
        rho_next[-1] = 2 * rho_next[-2] - rho_next[-3]
        V_next[-1] = 2 * V_next[-2] - V_next[-3]
        T_next[-1] = 2 * T_next[-2] - T_next[-3]
        
        
        # Update data containers
        rho = rho_next.copy()
        V = V_next.copy()
        T = T_next.copy()
        
        rho_n[t] = rho.copy()
        V_n[t] = V.copy()
        T_n[t] = T.copy()
        
        t_n[t] = dt
    
    return t_n, rho_n, V_n, T_n

def plot_flow_properties(x, p, rho, T, p_an, rho_an, T_an):
    pyplot.figure(figsize=(16,14))
    pyplot.subplot(311)
    pyplot.plot(x, p, color='k', ls='-', lw=3, label='Numerical')
    pyplot.plot(x, p_an, color='r', ls='--', lw=3, label='Analytical')
    pyplot.ylim(0,1)
    #pyplot.xlabel(r'$x$' + '\nNondimensional distance through nozzle')
    pyplot.ylabel('Nondimensional Pressure\n' + r'$\frac{p}{p_0}$')
    pyplot.legend();
    
    pyplot.subplot(312)
    pyplot.plot(x, rho, color='k', ls='-', lw=3, label='Numerical')
    pyplot.plot(x, rho_an, color='r', ls='--', lw=3, label='Analytical')
    pyplot.ylim(0,1)
    #pyplot.xlabel(r'$x$' + '\nNondimensional distance through nozzle')
    pyplot.ylabel('Nondimensional Density\n' + r'$\frac{\rho}{\rho_0}$')
    pyplot.legend();
    
    pyplot.subplot(313)
    pyplot.plot(x, T, color='k', ls='-', lw=3, label='Numerical')
    pyplot.plot(x, T_an, color='r', ls='--', lw=3, label='Analytical')
    pyplot.ylim(0,1)
    pyplot.xlabel(r'$x$' + '\nNondimensional distance through nozzle')
    pyplot.ylabel('Nondimensional Temperature\n' + r'$\frac{T}{T_0}$')
    pyplot.legend();
