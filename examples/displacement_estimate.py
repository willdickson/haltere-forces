import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from haltere_forces.halteres import Halteres
from parameters import drosophila_param
from parameters import calliphora_param

#param = calliphora_param
param = drosophila_param

for k, v in param.items():
    print(k,v)

def m_to_mm(val):
    return val*1e3

def m_to_um(val):
    return val*1e6

param['shift'] = 0.0
param['tilt_angle'] = 0.0

hsim = Halteres(param=param)

t = hsim.t
dt = t[1] - t[0]
angle = hsim.angle
omega_deg = np.array([1000.0,   0.0,   0.0])
omega = np.deg2rad(omega_deg)
force = hsim.force(omega)
coriolis = force['left']['lateral']['coriolis']

acc = -coriolis/drosophila_param['mass']
vel = sp.integrate.cumulative_trapezoid(acc,dx=dt,initial=0.0)
pos = sp.integrate.cumulative_trapezoid(vel,dx=dt, initial=0.0)

if 0:
    fig, ax = plt.subplots(3,1)
    ax[0].plot(t, m_to_um(acc))
    ax[0].set_ylabel(r'acc ($\mu m$/$s^2$)')
    ax[0].grid(True)
    
    ax[1].plot(t, m_to_um(vel))
    ax[1].set_ylabel(r'vel ($\mu m$/$s$)')
    ax[1].grid(True)
    
    ax[2].plot(t, m_to_um(pos))
    ax[2].set_ylabel(r'pos ($\mu m$)')
    ax[2].set_xlabel('t ($s$)')
    ax[2].grid(True)

if 1:
    fig, ax = plt.subplots(4,1)

    omega_str = np.array2string(omega_deg, precision=2, separator=',')
    ax[0].set_title(r'$\omega$ = ' + omega_str + '  (deg/s)')
    
    ax[0].plot(t, m_to_um(pos))
    ax[0].set_ylabel(r'body pos ($\mu m$)')
    ax[0].set_xlabel('t ($s$)')
    ax[0].grid(True)

    ax[1].plot(t, m_to_um(acc))
    ax[1].set_ylabel(r'body acc ($\mu m/s^2$)')
    ax[1].set_xlabel('t ($s$)')
    ax[1].grid(True)
    
    ax[2].plot(t, m_to_um(hsim.pos_left[:,2]))
    ax[2].set_ylabel(r'haltere vert pos ($\mu m$)')
    ax[2].grid(True)

    ax[3].plot(t, coriolis)
    ax[3].grid(True)
    ax[3].set_ylabel(r'coriolis ($N$)')
    ax[3].set_xlabel('t ($s$)')
plt.show()
