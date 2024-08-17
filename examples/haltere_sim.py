import numpy as np
import matplotlib.pyplot as plt
from haltere_forces.halteres import Halteres


param = {
        'num_pt'      : 10000,    
        'num_cycle'   : 3, 
        'waveform'    : 'filtered_triangle',
        'cutoff_freq' : 130.0,               # (Hz) 
        'frequency'   : 130.0,               # (Hz)
        'amplitude'   : np.deg2rad(90.0),    # (deg)
        'mass'        : 4.75e-9,             # (kg) 
        'length'      : 0.9e-3,              # (m)
        'separation'  : 2.5e-3,              # (m)
        'tilt_angle'  : np.deg2rad(30.0),    # (rad)
        }

h = Halteres(param=param)
pos_l, pos_r = h.positions

print(pos_l[:10,:])
print()
print(pos_r[:10,:])

fig, ax = plt.subplots(1,1)
ax.plot(h.t, h.angle)
ax.set_xlabel('t')
ax.set_ylabel('angle')
ax.grid(True)
plt.show()


