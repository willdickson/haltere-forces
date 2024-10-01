import numpy as np
import matplotlib.pyplot as plt
from haltere_forces import waveform

# Triangle waveform parameters
num_pt = 10000
num_cycle = 2
freq = 200.0
period = 1/freq 
amplitude = 90.0 
cutoff_freq = 2*freq 
rescale = True 

# Create unfiltered and filtered triangle waves
t = np.linspace(0, num_cycle*period, num_pt)
x = waveform.triangle(t, amplitude, period)
t, y = waveform.filtered_triangle(
        num_pt=num_pt,
        num_cycle=num_cycle, 
        amplitude=amplitude, 
        period=period, 
        cutoff_frequency=cutoff_freq, 
        rescale=rescale,
        )

# Compute derivatives
dt = t[1]-t[0]
dx = np.gradient(x)/dt
dy = np.gradient(y)/dt

# Plot waveforms
fig, ax = plt.subplots(2,1, sharex=True, figsize=(8,7))

line_flt, = ax[0].plot(t, y, 'b')
ax[0].set_ylabel(r'haltere angle, $\theta$,  (deg)')
ax[0].grid(True)
#ax[0].set_title('Haltere angle')
#ax[0].set_ylim(-amplitude, amplitude)

ax[1].plot(t, dy, 'b')
ax[1].set_ylabel(r' $\dot{\theta} (deg/s)$ ')
ax[1].grid(True)

ax[1].set_xlim(0, num_cycle*period)


#ax[1].set_xlabel('t (sec)')
fig.tight_layout()
fig.savefig('haltere_angles.png')
plt.show()
