#!/usr/bin/env python2.6

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import numpy as np
import sys, math



yrange = None
for arg in sys.argv:
  if arg.startswith('-y'):
    yrange = [float(x) for x in arg[2:].split(',')]


for filename in sys.argv[1:]:
  dta = 0.
  num = 0
  sid = -1
  x = []
  y = []
  for line in open(filename, 'r'):
    sp = line.split()
    if line.startswith('* Defining session'):
      x.append([])
      y.append([])
      sid += 1
    if len(sp) > 2 and sp[2] == 'event':
      try:
        sid = int(sp[0].strip('#')) - 1
        dwell_total = float(sp[-1].split('=')[1][:-3])
        if dwell_total > 1.0e9: continue
        dta += dwell_total
        num += 1
        x[sid].append(num)
        y[sid].append(dta / num)
      except ValueError:
        pass
      except IndexError:
        pass

  window = 500
  for i in range(len(x)):
    X = x[i]
    Y = y[i]
    label = 'Session %d Total - t = %.1e' % (i+1, Y[-1])
    plt.plot(X, Y, label=label)

plt.legend(loc='upper right', prop={'size':12})
plt.ylabel('Average Total Non-catalytic Association Time', fontsize=12)
plt.xlabel('Completed Substrate Replicate Simulations', fontsize=12)
if yrange != None:
  plt.ylim(yrange)

plt.title(sys.argv[1])
plt.show()
#fig = matplotlib.pyplot.gcf()
#fig.savefig(sys.argv[2], dpi=300)
