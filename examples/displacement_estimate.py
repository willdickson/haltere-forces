import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from haltere_forces.halteres import Halteres
from parameters import drosophila_param

def m_to_mm(val):
    return val*1e3

body_length = 3.0e-3
drosophila_param['shift'] = 0.0
drosophila_param['tilt_angle'] = 0.0

hsim = Halteres(param=drosophila_param)

t = hsim.t
dt = t[1] - t[0]
angle = hsim.angle
omega = np.deg2rad(np.array([1000.0,   0.0,   0.0]))
force = hsim.force(omega)
coriolis = force['left']['lateral']['coriolis']


acc = coriolis/drosophila_param['mass']
vel = sp.integrate.cumulative_trapezoid(acc,dx=dt,initial=0.0)
pos = sp.integrate.cumulative_trapezoid(vel,dx=dt, initial=0.0)

fig, ax = plt.subplots(3,1)
ax[0].plot(t, m_to_mm(acc))
ax[0].set_ylabel(r'acc ($mm$/$s^2$)')
ax[0].grid(True)

ax[1].plot(t, m_to_mm(vel))
ax[1].set_ylabel(r'vel ($mm$/$s$)')
ax[1].grid(True)

ax[2].plot(t, m_to_mm(pos))
ax[2].set_ylabel('pos ($mm$)')
ax[2].set_xlabel('t ($s$)')
ax[2].grid(True)
plt.show()
