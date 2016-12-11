import numpy
from matplotlib import pyplot
import math
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] =16
from matplotlib import animation
from IPython.display import HTML


nx =101
dx =0.01
nt= 150
x = numpy.linspace(-20,20,nx)
D =1
alpha = 1
u = 8*(numpy.exp(-5*x**2)) 
fig = pyplot.figure(figsize=(8,6))
ax = pyplot.axes(xlim=(-20,20),ylim=(-0.5,9))
ax.set_xlabel('$x$')
ax.set_ylabel('$u$')
#ax.text(0,-0.25,'See how the curve levels off at u =1',color='red')
ax.annotate('See how the curve levels off at u =1',color='#ff7f0e',xy=(-5,1),xytext=(-9,1.4), arrowprops =dict(facecolor='black', shrink =0.05))
line = ax.plot([],[],color='c',ls=':',lw=3)[0]


def CatchtheFisherwaveanimate(n):
    dt =0.0001
    un= numpy.zeros_like(u)
    for n in range (nt):
        un= numpy.copy(u)
        u[1:-1] = un[1:-1] \
                + D * dt/dx**2 * (un[2:] - 2*un[1:-1] + un[0:-2]) \
                + alpha * un[1:-1] * dt * (1 - un[1:-1])       
        u[0] = u[1] 
        u[-1] = u[-2]
        line.set_data(x,u)
anim = animation.FuncAnimation(fig,CatchtheFisherwaveanimate,\
                               frames=nt,interval=60)