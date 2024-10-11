import numpy as np
import matplotlib.pyplot as plt
from haltere_forces import waveform
from haltere_forces.simplified_coriolis import lateral_coriolis_from_pitch
from parameters import drosophila_param as param

# Triangle waveform parameters
num_pt = 10000
num_cycle = 2
freq = param['frequency']
period = 1/freq
amplitude = param['amplitude'] 
cutoff_freq = param['cutoff_freq'] 
omega_x = np.deg2rad(300.0)
rescale = True 

# Create unfiltered and filtered triangle waves
t = np.linspace(0, num_cycle*period, num_pt)
t, theta = waveform.filtered_triangle(
        num_pt=num_pt,
        num_cycle=num_cycle, 
        amplitude=amplitude, 
        period=period, 
        cutoff_frequency=cutoff_freq, 
        rescale=rescale,
        )

# Compute derivatives
dt = t[1]-t[0]
dtheta_dt = np.gradient(theta)/dt

f = lateral_coriolis_from_pitch(
        param['mass'], 
        param['length'], 
        param['tilt_angle'],
        omega_x, 
        theta, 
        dtheta_dt
        )

# Plot waveforms
fig, ax = plt.subplots(3,1, sharex=True, figsize=(8,7))

line_flt, = ax[0].plot(t, theta, 'b')
#ax[0].set_ylabel(r'$\theta$, (deg)')
ax[0].grid(True)

ax[1].plot(t, dtheta_dt, 'b')
#ax[1].set_ylabel(r'$\dot{\theta} (deg/s)$ ')
ax[1].grid(True)

ax[2].plot(t, f, 'b')
#ax[2].set_ylabel('lateral coriolis (N)')
ax[2].grid(True)

ax[2].set_xlim(0, num_cycle*period)
#ax[2].set_xlabel('t (sec)')

fig.tight_layout()
fig.savefig('lateral_coriolis_forces.png')
plt.show()
