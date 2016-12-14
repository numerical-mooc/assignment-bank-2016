import numpy
from matplotlib import pyplot
import math
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] =16
from matplotlib import animation
from IPython.display import HTML


nx =201
dx =0.2
nt= 300
x = numpy.linspace(0,10,nx)
D =1
alpha = 1

u= numpy.zeros(nx)
lbound = numpy.where(x <= 1)
ubound = numpy.where(x >= 0)            
bounds = numpy.intersect1d(lbound, ubound)
u[bounds]=1       

fig = pyplot.figure(figsize=(8,6))
ax = pyplot.axes(xlim=(0,10),ylim=(0,1.5))
ax.set_xlabel('$x$')
ax.set_ylabel('$u$')
#ax.text(0,-0.25,'See how the curve levels off at u =1',color='red')
ax.annotate('See how the curve',color='black',wrap= True, xy=(7,1),xytext=(3,1.3))
ax.annotate('levels off at u =1',color='black',wrap= True,xy=(7,1),xytext=(3,1.2),arrowprops =dict(facecolor='black', shrink =0.04))
line = ax.plot([],[],color='c',ls=':',lw=3)[0]



def CatchtheFisherwaveanimate(n):
    dt =0.0003
    un= numpy.zeros_like(u)
    ax.text(6,0.3,'Absence of gene',color='green')
    ax.fill_between(x,un,1.5, color='aqua')
    ax.fill_between(x,u,0, color='green')
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