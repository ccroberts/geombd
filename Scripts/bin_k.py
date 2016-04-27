#!/usr/bin/python

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys

k = []

filename = sys.argv[1]
binSize = float(sys.argv[2])
serach = sys.argv[3]
print serach
for j,line in enumerate(open(filename, 'r')):
  sp = line.split()
  if line.strip().startswith(serach):
    if sp[2] == 'k_on':
      k.append(float(sp[4]))
    if sp[4] == 'k_on':
      k.append(float(sp[6]))

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
