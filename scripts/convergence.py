#!/usr/bin/env python2.6

import matplotlib
matplotlib.use('Agg')
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

print sys.argv[1]

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
    #bta = bnd/num
    x[sid][0].append(num)
    y[sid][0].append(bnd/num)
  if line.startswith("   (session") and sp[2] == 'bs':
    sid = int(sp[1]) - 1
    bsid = int(sp[3][:-1]) + 1
    kon = float(sp[6])
    bnd = float(sp[8].split('=')[1])
    num = int(sp[9].split('=')[1])
    #bta = bnd/num
    x[sid][bsid].append(num)
    y[sid][bsid].append(bnd/num)

def running_variance(x, y, w):
  X = []
  Y = []
  for i in range(0, len(y)-w):
    yi = y[i:(i+w)]
    m = sum(yi)/len(yi)
    v = math.sqrt(sum([pow(yii-m,2) for yii in yi]))/len(yi)
    X.append(x[i+w-1])
    Y.append(v)
  return X, Y


window = 200
if b_stotal:
  for ci in range(len(r_stotal) - 1):
    newy = []
    for i in range(len(x[0][0])):
      newyi = 0.
      for j in range(r_stotal[ci], r_stotal[ci+1]):
        newyi += y[0][j+1][i]
      newy.append(newyi)
    if len(newy) < window:
      window = len(newy)
    XX, YY = running_variance(x[0][0], newy, window)
    label = 'Combined Session %d Total - c = %.1e' % (ci+1, sum(YY[-window:])/len(YY[-window:]))
    plt.plot(XX, YY, label=label)
else:
  for i in range(len(x)):
    for j in range(len(x[i])):
      if len(y[i][j]) < window:
        window = len(y[i][j]) - 1
      if j == 0 and b_total:
        X, Y = running_variance(x[i][j], y[i][j], window)
        label = 'Session %d Total - c = %.1e' % (i+1, sum(Y[-window:])/len(Y[-window:]))
        plt.plot(X, Y, label=label)
      if j > 0 and b_bs:
        X, Y = running_variance(x[i][j], y[i][j], window)
        label = 'Session %d BS %d - c = %.1e' % (i+1, j-1, sum(Y[-window:])/len(Y[-window:]))
        plt.plot(X, Y, label=label)

if yrange != None:
  plt.ylim(yrange)
plt.legend(loc='upper right', prop={'size':12})
plt.ylabel('Bound Fraction Variance', fontsize=16)
plt.xlabel('Completed Substrate Replicate Simulations', fontsize=16)

plt.title(sys.argv[1])
#plt.show()
fig = matplotlib.pyplot.gcf()
fig.savefig(sys.argv[2], dpi=300)
