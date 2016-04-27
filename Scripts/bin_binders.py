#!/usr/bin/python

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys

k = []

filename = sys.argv[1]
binSize = float(sys.argv[2])

for j,line in enumerate(open(filename, 'r')):
  sp = line.split()
  if len(sp) > 2 and sp[1] == 'Binding' and sp[2] == 'event':
    try:
      tb = float(sp[4].split('=')[1])
      k.append(tb / 1000.) #ps to ns
    except IndexError:
      pass


nbin = int((max(k)-min(k))/binSize)
hist, bins = np.histogram(k, bins=nbin)
width = 0.65 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2

plt.bar(center, hist, align='center', width=width)

#plt.legend(loc='upper left', prop={'size':9})
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
#plt.ylim([0, 100.])
#plt.tight_layout();
plt.show()

#fig = matplotlib.pyplot.gcf()
#fig.set_size_inches(16,12)
#fig.savefig("figure.png", dpi=300)
