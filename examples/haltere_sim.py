import numpy as np
import matplotlib.pyplot as plt
from haltere_forces.halteres import Halteres
from haltere_forces.halteres import reshape_to_nx3


#param = {
#        'num_pt'      : 100,    
#        'num_cycle'   : 1, 
#        'waveform'    : 'filtered_triangle',
#        'cutoff_freq' : 2*130.0,             # (Hz) 
#        'frequency'   : 130.0,               # (Hz)
#        'amplitude'   : np.deg2rad(90.0),    # (deg)
#        'mass'        : 5.89e-9,             # (kg) 
#        'length'      : 0.9e-3,              # (m)
#        'separation'  : 2.5e-3,              # (m)
#        'tilt_angle'  : np.deg2rad(30.0),    # (rad)
#        }

param = {
        'num_pt'      : 100,    
        'num_cycle'   : 1, 
        'waveform'    : 'filtered_triangle',
        'cutoff_freq' : 2*130.0,             # (Hz) 
        'frequency'   : 130.0,               # (Hz)
        'amplitude'   : np.deg2rad(90.0),    # (deg)
        'mass'        : 5.89e-9,             # (kg) 
        'length'      : 0.9e-3,              # (m)
        'separation'  : 2.5e-3,              # (m)
        'tilt_angle'  : np.deg2rad(30.0),    # (rad)
        }

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
    ax[0].set_title(r'$\omega$ = ' + omega_str)
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


