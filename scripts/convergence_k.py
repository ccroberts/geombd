#!/usr/bin/env python2.6

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import numpy as np
import sys, math



b_total = True
b_bs = True
yrange = None
for arg in sys.argv:
  if arg == '-nototal':
    b_total = False
  elif arg == '-total':
    b_bs = False
  elif arg.startswith('-y'):
    yrange = [float(x) for x in arg[2:].split(',')]


sid = -1
x = []
y = []
for line in open(sys.argv[1], 'r'):
  sp = line.split()
  if line.startswith('* Defining session'):
    x.append([ [] ])
    y.append([ [] ])
    sid += 1
  if line.startswith(' + Binding criteria'):
    x[sid].append([])
    y[sid].append([])
  if line.startswith("   (session") and sp[1][-1] == ')':
    sid = int(sp[1][:-1]) - 1
    kon = float(sp[4])
    bnd = float(sp[6].split('=')[1])
    num = float(sp[7].split('=')[1])
    bta = bnd/num
    x[sid][0].append(num)
    y[sid][0].append(kon)
  if line.startswith("   (session") and sp[2] == 'bs':
    sid = int(sp[1]) - 1
    bsid = int(sp[3][:-1]) + 1
    kon = float(sp[6])
    bnd = float(sp[8].split('=')[1])
    num = int(sp[9].split('=')[1])
    bta = bnd/num
    x[sid][bsid].append(num)
    y[sid][bsid].append(kon)

for i in range(len(x)):
  for j in range(len(x[i])):
    if j == 0 and b_total:
      X = x[i][j]
      Y = y[i][j]
      label = 'Session %d Total - k = %.1e' % (i+1, sum(Y[-500:])/500.)
      plt.plot(X, Y, label=label)
    if j > 0 and b_bs:
      X = x[i][j]
      Y = y[i][j]
      label = 'Session %d BS %d - k = %.1e' % (i+1, j-1, sum(Y[-500:])/500.)
      plt.plot(X, Y, label=label)

plt.legend(loc='upper right', prop={'size':12})
plt.ylabel('Rate Constant Variance', fontsize=16)
plt.xlabel('Completed Substrate Replicate Simulations', fontsize=16)
if yrange != None:
  plt.ylim(yrange)

plt.title(sys.argv[1])
#plt.show()
fig = matplotlib.pyplot.gcf()
fig.savefig(sys.argv[2], dpi=300)
