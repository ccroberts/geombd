#!/usr/bin/env python2.6

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import numpy as np
import sys


x = []
y = []


for line in open(sys.argv[1], 'r'):
  if line.startswith("   (session "):
    sp = line.split()
    sid = int(sp[1][:-1]) - 1
    kon = float(sp[4])
    num = int(sp[7].split('=')[1])
    try:
      x[sid].append(num)
      y[sid].append(kon)
    except IndexError:
      x.append([])
      y.append([])
      x[sid].append(num)
      y[sid].append(kon)

maxy = 0.
for i in range(len(x)):
  plt.plot(x[i], y[i])
  imaxy = max(y[i])
  maxy = max([maxy, imaxy])

#plt.legend(loc='upper left', prop={'size':9})
#plt.xlabel('Association Time (ns, %dns bins)' % binSize, fontsize=18)
#plt.ylabel('Percent Ligands Associated (%)', fontsize=18)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
#plt.tight_layout();
plt.show()

#fig = matplotlib.pyplot.gcf()
#fig.set_size_inches(16,12)
#fig.savefig(sys.argv[2], dpi=300)
