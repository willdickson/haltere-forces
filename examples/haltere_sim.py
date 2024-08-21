import numpy as np
import matplotlib.pyplot as plt
from haltere_forces.halteres import Halteres
from haltere_forces.halteres import reshape_to_nx3


param = {
        'num_pt'      : 1000,    
        'num_cycle'   : 1, 
        'waveform'    : 'filtered_triangle',
        'cutoff_freq' : 2*130.0,             # (Hz) 
        'frequency'   : 130.0,               # (Hz)
        'amplitude'   : np.deg2rad(90.0),    # (deg)
        'mass'        : 5.89e-9,             # (kg) 
        'length'      : 0.9e-3,              # (m)
        'separation'  : 2.5e-3,              # (m)
        'tilt_angle'  : np.deg2rad(0.0),     # (rad)
        }

hsim = Halteres(param=param)

omega = np.array([0.0, 320.0, 0.0])
omega = np.deg2rad(omega)

t = hsim.t
ang = hsim.angle
pos = hsim.kinematics['left']['pos']
pos_norm = np.linalg.norm(pos,axis=1)
pos_norm = pos_norm[..., np.newaxis].repeat(3,1)
pos_unit = pos/pos_norm

force = hsim.force(omega)
coriolis = force['left']['coriolis']
angular_acc = force['left']['angular_acc']

coriolis_radial = np.sum(coriolis*pos_unit,axis=1)


fig, ax = plt.subplots(2,1,sharex=True)
ax[0].plot(t, ang)
ax[0].set_ylabel('angle (rad)')
ax[0].grid(True)

#ax[1].plot(t, coriolis[:,2])
ax[1].plot(t, coriolis_radial)
ax[1].set_ylabel('f_coriolis (N)')
ax[1].grid(True)

ax[1].set_xlabel('t (s)')
plt.show()


