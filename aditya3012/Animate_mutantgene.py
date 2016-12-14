import numpy
from matplotlib import pyplot
import math
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] =16
from matplotlib import animation
from IPython.display import HTML


nx =200
dx =0.1
nt= 350
x = numpy.linspace(-10,10,nx)
D =1
alpha = 1

u = numpy.zeros(nx)
lbound = numpy.where(x >= -2)
ubound = numpy.where(x <= 2)
bounds = numpy.intersect1d(lbound, ubound)
u[bounds]=0.01

fig = pyplot.figure(figsize=(8,6))
ax = pyplot.axes(xlim=(-10,10),ylim=(0,1.5))
ax.set_xlabel('$x$')
ax.set_ylabel('$u$')
#ax.text(0,-0.25,'See how the curve levels off at u =1',color='red')
line = ax.plot([],[],color='c',ls=':',lw=3)[0]

def CatchtheFisherwaveanimate(n):
    dt =0.0001
    un= numpy.zeros_like(u)
    ax.fill_between(x,u,0,interpolate=True, color='blue')
    ax.text(-2.5,.8,'Sickle cell gene',color='white')
    for n in range (nt):
        un= numpy.copy(u)
        u[1:-1] = un[1:-1] \
                + D * dt/dx**2 * (un[2:] - 2*un[1:-1] + un[0:-2]) \
                + alpha * un[1:-1] * dt * (1 - un[1:-1])       
        u[0] = u[1] 
        u[-1] = u[-2]
        line.set_data(x,u)
anim = animation.FuncAnimation(fig,CatchtheFisherwaveanimate,\
                               frames=nt,interval=50)