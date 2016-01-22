#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys

k = []

filename = sys.argv[1]
binSize = float(sys.argv[2])
fileout = sys.argv[3]

for j,line in enumerate(open(filename, 'r')):
  sp = line.split()
  if len(sp) > 2 and sp[1] == 'Binding' and sp[2] == 'event':
    try:
      tl = float(sp[6].split('=')[1][:-3])
      tm = float(sp[7].split('=')[1][:-3])
      tt = float(sp[7].split('=')[1][:-3])
      k.append(tt) #ps to ns
    except IndexError:
      pass
    except ValueError:
      pass


nbin = int(max(k)/binSize) + 1
hist, bins = np.histogram(k, bins=nbin)
total = float(sum(hist))
hist = [100.*float(x)/total for x in hist]
width = 0.9*binSize
center = list(range(0, int(nbin*binSize), binSize))
center = [x+(binSize/2.0) for x in center]

print len(hist), len(center)
plt.bar(center, hist, align='center', width=width)

plt.ylim([0, 50.])
plt.yticks(fontsize=12)
plt.gca().yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f%%'))

plt.xlim([0, 25.0e3])
plt.xticks(fontsize=12)
plt.axes().ticklabel_format(style='sci', axis='x', scilimits=(0,0))

plt.show()

fig = matplotlib.pyplot.gcf()
#fig.set_size_inches(16,12)
fig.savefig(fileout, dpi=300)
