#!/usr/bin/env python2.6

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import numpy as np
import sys, math



b_total = True
b_bs = True
b_stotal = False
r_stotal = []
yrange = None
for arg in sys.argv:
  if arg == '-nototal':
    b_total = False
  elif arg == '-total':
    b_bs = False
  elif arg.startswith('-c'):
    b_stotal = True
    r_stotal = [int(x) for x in arg[2:].split(',')]
    r_stotal.insert(0, 0)
  elif arg.startswith('-y'):
    yrange = [float(x) for x in arg[2:].split(',')]


sid = -1
bnd = 0
x = []
y = []
b = []
for line in open(sys.argv[1], 'r'):
  sp = line.split()
  if line.startswith('* Defining session'):
    x.append([ [] ])
    y.append([ [] ])
    b.append([ [] ])
    sid += 1
  if line.startswith(' + Binding criteria'):
    x[sid].append([])
    y[sid].append([])
    b[sid].append([])
  if line.startswith("   (session") and sp[1][-1] == ')':
    sid = int(sp[1][:-1]) - 1
    kon = float(sp[4])
    bnd = float(sp[6].split('=')[1])
    num = float(sp[7].split('=')[1])
    bta = bnd/num
    x[sid][0].append(num)
    y[sid][0].append(kon)
    b[sid][0].append(bnd)
  if line.startswith("   (session") and sp[2] == 'bs':
    sid = int(sp[1]) - 1
    bsid = int(sp[3][:-1]) + 1
    kon = float(sp[6])
    bnd = float(sp[8].split('=')[1])
    num = int(sp[9].split('=')[1])
    bta = bnd/num
    x[sid][bsid].append(num)
    y[sid][bsid].append(kon)
    b[sid][bsid].append(bnd)

window = 100
if b_stotal:
  for ci in range(len(r_stotal) - 1):
    newy = []
    for i in range(len(x[0][0])):
      newyi = 0.
      inewyi = 0
      for j in range(r_stotal[ci], r_stotal[ci+1]):
        newyi += y[0][j+1][i]
        inewyi += 1
      #newyi /= inewyi
      newy.append(newyi)
    label = 'Combined Session %d Total - k = %.1e +/- %.1f%%' % (ci+1, sum(newy[-window:])/len(newy[-window:]), pow(b[0][-1], -0.5) * 100.)
    plt.plot(x[0][0], newy, label=label)
else:
  for i in range(len(x)):
    for j in range(len(x[i])):
      if j == 0 and b_total:
        X = x[i][j]
        Y = y[i][j]
        label = 'Session %d Total - k = %.1e +/- %.1f%%' % (i+1, sum(Y[-window:])/len(Y[-window:]), pow(b[i][j][-1], -0.5)*100.)
        plt.plot(X, Y, label=label)
      if j > 0 and b_bs:
        X = x[i][j]
        Y = y[i][j]
        label = 'Session %d BS %d - k = %.1e +/- %.1f%%' % (i+1, j-1, sum(Y[-window:])/len(Y[-window:]), pow(b[i][j][-1], -0.5)*100.)
        plt.plot(X, Y, label=label)
        for i in range(0, len(X), 2):
          print X[i], Y[i]

plt.legend(loc='upper right', prop={'size':12})
plt.ylabel('Rate Constant', fontsize=16)
plt.xlabel('Completed Substrate Replicate Simulations', fontsize=16)
if yrange != None:
  plt.ylim(yrange)

plt.title(sys.argv[1])
plt.show()
#fig = matplotlib.pyplot.gcf()
#fig.savefig(sys.argv[2], dpi=300)
