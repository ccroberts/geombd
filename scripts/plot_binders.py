#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys

def usage():
  print 'Not enough arguments! Usage:'
  print '\t%s BINSIZE_ns LOGFILE1...' % sys.argv[0]
  sys.exit()

if len(sys.argv) < 3:
  usage()

try:
  binSize = int(sys.argv[1])
except ValueError:
  usage()
  

T = []
N = []

filename = sys.argv[2]
for j,line in enumerate(open(filename, 'r')):
  sp = line.split()
  if line.startswith(' + Ligand replicates:'):
    N.append(int(sp[-1]))
    T.append([])
  try:
    if sp[1] == 'Binding' and sp[2] == 'event':
      n = int(sp[0][1:]) - 1
      t = float(sp[4].split('=')[1]) / 1e3 # ps -> ns
      T[n].append(t)
  except IndexError:
    pass
  except ValueError:
    pass

for i in range(len(N)):
  # bin the binding times
  maxt = max(T[i])
  nbin = int(round(maxt) / binSize)
  hist, bins = np.histogram(T[i], bins=nbin)

  hist = [100.*(float(x)/float(N[i])) for x in list(hist)]
  width = 0.65 * (bins[1] - bins[0])
  center = (bins[:-1] + bins[1:]) / 2

  plt.bar(center, hist, align='center', width=width)

'''
  # cumulative
  x = []
  y = []

  tt = 0.
  for i,h in enumerate(hist):
    tt += h
    y.append(tt)
    x.append(bins[i])

  plt.plot(x, y)
'''

#plt.legend(loc='upper left', prop={'size':9})
plt.xlabel('Association Time (ns, %dns bins)' % binSize, fontsize=18)
plt.ylabel('Percent Ligands Associated (%)', fontsize=18)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
#plt.ylim([0, 100.])
#plt.tight_layout();
#plt.show()

fig = matplotlib.pyplot.gcf()
#fig.set_size_inches(16,12)
fig.savefig("figure.png", dpi=300)
