#!/usr/bin/env python2.6

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
import numpy as np
import sys, math



yrange = None
for arg in sys.argv:
  if arg.startswith('-y'):
    yrange = [float(x) for x in arg[2:].split(',')]


for filename in sys.argv[1:]:
  dta = 0.
  num = 0
  sid = -1
  x = []
  y = []
  A = 0
  B = []
  for line in open(filename, 'r'):
    sp = line.split()
    if line.startswith('* Defining session'):
      x.append([])
      y.append([])
      sid += 1
    if len(sp) > 2 and sp[2] == 'event':
      if line.find('event', line.find('event')+1) > -1:
        continue
      try:
        sid = int(sp[0].strip('#')) - 1
        #dwell_temp = float(sp[-3].split('=')[1][:-3])
        dwell_total = float(sp[-2].split('=')[1][:-3])
        dta += dwell_total
        num += 1
        x[sid].append(num)
        y[sid].append(dta / num)
        if dwell_total != 0.: A += 1
        B.append(dwell_total)
      except ValueError:
        pass
      except IndexError:
        pass

  window = 500
  for i in range(len(x)):
    X = x[i]
    Y = y[i]
    print str('%30s' % filename), ':: avg_total=', Y[-1], '\tavg_nca=', sum(B)/A, '\tmax =', max(B), '\tstdev=', np.std(B), '\tstderr=', (np.std(B) / math.sqrt(A)),'\tNnca=', A, '\tN=', num, '\t%nca=', 100. * float(A)/float(num)
    #label = 'Session %d Total - t = %.1e +/- %.1f%%' % (i+1, Y[-1], pow(A, -0.5) * 100.)
    #plt.plot(X, Y)#, label=label)

'''
plt.legend(loc='upper right', prop={'size':12})
plt.ylabel('Average Total Non-catalytic Association Time', fontsize=12)
plt.xlabel('Completed Substrate Replicate Simulations', fontsize=12)
if yrange != None:
  plt.ylim(yrange)

plt.title(sys.argv[1])
#fig = matplotlib.pyplot.gcf()
#fig.savefig(sys.argv[2], dpi=300)
'''
