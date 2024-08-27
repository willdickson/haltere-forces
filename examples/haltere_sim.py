import numpy as np
import matplotlib.pyplot as plt
from haltere_forces.halteres import Halteres
from haltere_forces.halteres import reshape_to_nx3
from parameters import drosophila_param
from parameters import calliphora_param

param = drosophila_param

for k,v in param.items():
    print(k,v)


hsim = Halteres(param=param)

t = hsim.t
ang = hsim.angle

omega_list = [
        np.array([320.0,   0.0,   0.0]),
        np.array([  0.0, 320.0,   0.0]),
        np.array([  0.0,   0.0, 320.0]),
        ]

for omega_deg in omega_list:

    omega = np.deg2rad(omega_deg)
    force = hsim.force(omega)

    fig, ax = plt.subplots(3,1,sharex=True)
    omega_str = np.array2string(omega_deg, precision=2, separator=',')
    ax[0].set_title(r'$\omega$ = ' + omega_str + '  (deg/s)')
    ax[0].plot(t, ang, 'b')
    ax[0].set_ylabel('angle (rad)')
    ax[0].grid(True)

    ax[1].plot(t, force['left']['lateral']['coriolis'], 'b')
    #ax[1].plot(t, force['left']['lateral']['total'], 'r')
    ax[1].set_ylabel('force left (N)')
    ax[1].grid(True)
    
    ax[2].plot(t, force['right']['lateral']['coriolis'], 'b')
    #ax[2].plot(t, force['right']['lateral']['total'], 'r')
    ax[2].set_ylabel('force right (N)')
    ax[2].grid(True)
    ax[2].set_xlabel('t (s)')

plt.show()


