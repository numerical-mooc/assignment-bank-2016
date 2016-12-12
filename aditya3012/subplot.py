import numpy
from matplotlib import pyplot
import math
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] =16


fig = pyplot.figure(figsize=(12,7))
fig.subplots_adjust(wspace=.3)
sub1 = fig.add_subplot(221)
sub2 = fig.add_subplot(222)
sub1.tick_params(labelsize=16),sub2.tick_params(labelsize=16)
sub2.set_ylim(-0.5,1.5)
sub1.set_xlabel('${x}$'), sub1.set_ylabel('${f(x) = x(1-x)}$') 
sub2.set_xlabel('${x}$'), sub2.set_ylabel('${u}$') 
