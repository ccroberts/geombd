#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot_sys(kdirect, Bdirect, kindirect, Bindirect, konfree, E1, label, alpha):
    alpha = 1.0
    kcat1 = 251.
    km1 = 0.026
    #E1 = 1e-9
    kindirect *= alpha
    konfree *= alpha
    S1 = 1e-3
    #kdirect = 1.7e5
    #Bdirect = 0.0011
    #kindirect = 6.9e7
    #Bindirect = 0.0002

    #konfree = 6.0e9

    kcat2 = 200.
    km2 = 0.4e-6
    E2 = E1

    r = []
    e = []

    for ps12 in [1e-9, 5e-9, 1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6,  1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 0.5, 1]:
      s12 = ps12
      vd = kdirect * E2
      vi = (1 - Bdirect) * kindirect * E2 * s12

      vf = konfree * E2 * s12

      r.append(s12)
      e.append((vd+vi)/vf)

    #plt.semilogx(r, e, label=label)
    plt.loglog(r, e, label=label)

'''
    for tp in range(1, 60*60*24):
      vd = kdirect * E1
      i_colc.append(i_colc[-1] + dt * (1 - Bdirect - Bindirect) * kcat1 * E1 * (s[-1] / (km1 + s[-1])))
      vi = (1 - Bdirect) * kindirect * E2 * i_colc[-1]
      maxv = kcat2 * E2 * (i_colc[-1] / (km2 + i_colc[-1]))
      #if vi+vd > maxv:
      #  p_colc.append(p_colc[-1] + (dt * maxv))
      #else:
      #  p_colc.append(p_colc[-1] + (dt * (vd + vi)))
      p_colc.append(p_colc[-1] + (dt * (vd + vi)))
      #i_colc[-1] = i_colc[-1] - (dt * (vd + vi))

      i_free.append(i_free[-1] + dt * kcat1 * E1 * (s[-1] / (km1 + s[-1])))
      vf = konfree * E2 * i_free[-1]
      maxv = kcat2 * E2 * (i_free[-1] / (km2 + i_free[-1]))
      #if vf > maxv: vf = maxv
      p_free.append(p_free[-1] + (dt * vf))
      #i_free[-1] = i_free[-1] - (dt * vf)

      t.append(tp/3600.)
      e.append(i_colc[-1] / i_free[-1])

    print label, e[-1]
    plt.plot(t, e, label=label)
'''

alpha = 1e-4
plot_sys(1.6e6, 0.0299, 6.3e9, 0.0228, 6.0e9, 1e-9, 'A.C. - 10nm', alpha)
plot_sys(2.2e6, 0.1431, 2.0e10, 0.0191, 6.0e9, 1e-9, 'Planar Origami - 10nm', alpha)
plot_sys(1.6e6, 0.1362, 1.9e10, 0.0180, 6.0e9, 1e-9, 'Tubular Origami - 10nm', alpha)

plt.legend(loc='upper right', prop={'size':12})
plt.xlabel("Concentration (M)", fontsize=16)
plt.ylabel("Colocalized Enhancement", fontsize=16)
#plt.ylim([0, 10000.])
#plt.show()
fig = matplotlib.pyplot.gcf()
fig.savefig('compare.png', dpi=300)
