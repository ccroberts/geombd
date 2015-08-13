#!/usr/bin/python

import sys, math

coords = []
center = [0., 0., 0.]
N = 0

for line in open(sys.argv[1], 'r'):
  if line.startswith('ATOM'):
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    coords.append([x,y,z])
    center[0] += x
    center[1] += y
    center[2] += z
    N += 1

center[0] /= N
center[1] /= N
center[2] /= N

rmsd = 0.

for c in coords:
  for i in range(3):
    dr = c[i] - center[i]
    rmsd += dr * dr

rmsd /= N
rmsd = math.sqrt(rmsd)
print('ROJ:', rmsd)
print('Center:', center)

