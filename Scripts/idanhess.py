#!/usr/bin/env python
import matplotlib
import matplotlib.pyplot as plt

def plot_sys(kdirect, Bdirect, kindirect, Bindirect, konfree, E1, label):
    kcat1 = 300.
    km1 = 0.026
    #E1 = 1e-9
    scale = 1.
    kindirect *= scale
    konfree *= scale
    S1 = 1e-3
    #kdirect = 1.7e5
    #Bdirect = 0.0011
    #kindirect = 6.9e7
    #Bindirect = 0.0002

    #konfree = 6.0e9

    kcat2 = 200.e20
    km2 = 0.4e-6
    E2 = E1

    t = [0]
    s = [S1]
    i_colc = [0]
    i_free = [0]
    p_colc = [0]
    p_free = [0]
    e = [0]

    dt = 1.0

    for tp in range(1, 60*60*24):
      vd = kdirect * E1
      i_colc.append(i_colc[-1] + dt * (1 - Bdirect - Bindirect) * kcat1 * E1 * (s[-1] / (km1 + s[-1])))
      vi = (1 - Bdirect) * kindirect * E2 * i_colc[-1]
      #maxv = kcat2 * E2 * (i_colc[-1] / (km2 + i_colc[-1]))
      #if vi+vd > maxv:
      #  p_colc.append(dt * maxv)
      #else:
      #  p_colc.append(dt * (vd + vi))
      p_colc.append(dt * (vd + vi))

      i_free.append(i_free[-1] + dt * kcat1 * E1 * (s[-1] / (km1 + s[-1])))
      vf = konfree * E2 * i_free[-1]
      #maxv = kcat2 * E2 * (i_free[-1] / (km2 + i_free[-1]))
      #if vf > maxv: vf = maxv
      p_free.append(dt * vf)

      t.append(tp/3600.)
      e.append(p_colc[-1] / p_free[-1])

    print label, e[-1]
    plt.plot(t, e, label=label)

plot_sys(1.6e6, 0.0299, 6.3e9, 0.0228, 6.0e9, 1e-7, 'A.C. - 10nm')
plot_sys(2.2e6, 0.1431, 2.0e10, 0.0191, 6.0e9, 1e-7, 'Planar Origami - 10nm')
plot_sys(1.6e6, 0.1362, 1.9e10, 0.0180, 6.0e9, 1e-7, 'Tubular Origami - 10nm')

plt.legend(loc='upper right', prop={'size':12})
plt.xlabel("Time (Hours)", fontsize=16)
plt.ylabel("Colocalized Enhancement", fontsize=16)
plt.ylim([0, 10.])
plt.show()
