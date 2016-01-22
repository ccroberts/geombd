#!/usr/bin/env python
import matplotlib
import matplotlib.pyplot as plt

def plot_sys(kdirect, Bdirect, kindirect, Bindirect, label):
    kcat1 = 300.
    km1 = 0.026
    E1 = 1e-9
    S1 = 1e-3
    #kdirect = 1.7e5
    #Bdirect = 0.0011
    #kindirect = 6.9e7
    #Bindirect = 0.0002

    konfree = 7.8e7

    kcat2 = 200.
    E2 = E1

    t = [0]
    s = [S1]
    i_colc = [0]
    i_free = [0]
    p_colc = [0]
    p_free = [0]
    e = [0]

    dt = 1.0

    for tp in range(1, 60*60*50):
      vd = kdirect * E1
      i_colc.append(i_colc[-1] + dt * (1 - Bdirect - Bindirect) * kcat1 * E1 * (s[-1] / (km1 + s[-1])))
      vi = (1 - Bdirect) * kindirect * E2 * i_colc[-1]
      if vi > kcat2: vi = kcat2
      p_colc.append(dt * (vd+vi))

      i_free.append(i_free[-1] + dt * kcat1 * E1 * (s[-1] / (km1 + s[-1])))
      vf = konfree * E2 * i_free[-1]
      if vf > kcat2: vf = kcat2
      p_free.append(dt * vf)

      t.append(tp/3600.)
      e.append(p_colc[-1] / p_free[-1])

    print label, e[-1]
    plt.plot(t, e, label=label)

plot_sys(1.7e5, 0.0011, 7.8e7, 0.0002, '10nm')
plot_sys(6.0e3, 0.0002, 7.8e7, 0.0002, '20nm')
plot_sys(1.5e2, 0.00002, 7.8e7, 0.0002, '45nm')
plot_sys(9.1e1, 0.000002, 7.8e7, 0.0002, '65nm')

plt.legend(loc='upper right', prop={'size':12})
plt.xlabel("Time (Hours)", fontsize=16)
plt.ylabel("Colocalized Enhancement", fontsize=16)
plt.ylim([0, 10.])
plt.show()
