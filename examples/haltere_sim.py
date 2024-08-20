import numpy as np
import matplotlib.pyplot as plt
from haltere_forces.halteres import Halteres
from haltere_forces.halteres import reshape_to_nx3


param = {
        'num_pt'      : 10,    
        'num_cycle'   : 1, 
        'waveform'    : 'filtered_triangle',
        'cutoff_freq' : 2*130.0,             # (Hz) 
        'frequency'   : 130.0,               # (Hz)
        'amplitude'   : np.deg2rad(90.0),    # (deg)
        'mass'        : 4.75e-9,             # (kg) 
        'length'      : 0.9e-3,              # (m)
        'separation'  : 2.5e-3,              # (m)
        'tilt_angle'  : np.deg2rad(30.0),    # (rad)
        }

h = Halteres(param=param)



#k = h.kinematics
#
#print(k['left']['pos'][:10,:])
#print()
#print(k['right']['pos'][:10,:])
#print(k['left']['pos'].shape)

omega = np.array([0, 10, 0])

if 1:
    fig, ax = plt.subplots(2,1)
    ax[0].plot(h.t, h.angle)
    ax[0].set_ylabel('angle')
    ax[0].grid(True)
    
    ax[1].plot(h.t, h.force_left(omega)[:,1])
    ax[1].set_ylabel('f')
    ax[1].set_xlabel('t')
    ax[1].grid(True)
    plt.show()


