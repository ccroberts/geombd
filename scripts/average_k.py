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
xstart = 0 
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
  elif arg.startswith('-x'):
    xstart = float(arg[2:])


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
    y4avg = []
    for i in range(len(newy)):
      if x[0][0][i] >= xstart:
        y4avg = newy[i:]
        break
    print 'Combined Session %d k = %.1e' % (ci+1, sum(y4avg)/len(y4avg))
else:
  for i in range(len(x)):
    for j in range(len(x[i])):
      if j == 0 and b_total:
        X = x[i][j]
        Y = y[i][j]
        y4avg = []
        for newi in range(len(Y)):
          if x[0][0][newi] >= xstart:
            y4avg = Y[newi:]
            break
        print 'Combined Session %d k = %.1e' % (i+1, sum(y4avg)/len(y4avg))
      if j > 0 and b_bs:
        X = x[i][j]
        Y = y[i][j]
        y4avg = []
        for newi in range(len(Y)):
          if x[0][0][newi] >= xstart:
            y4avg = Y[newi:]
            break
        print 'Combined Session %d k = %.1e' % (i+1, sum(y4avg)/len(y4avg))

