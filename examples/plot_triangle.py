import numpy as np
import matplotlib.pyplot as plt
from haltere_forces import waveform

# Triangle waveform parameters
num_pt = 10000
num_cycle = 3 
freq = 130.0
period = 1/freq 
amplitude = 5.0
cutoff_freq = freq 
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
fig, ax = plt.subplots(2,1, sharex=True)

line_tri, = ax[0].plot(t, x, 'b')
line_flt, = ax[0].plot(t, y, 'g')
ax[0].set_ylabel('x')
ax[0].grid(True)
ax[0].legend((line_tri, line_flt), ('no filter', 'filtered'), loc='upper right')
ax[0].set_title('Unfiltered and filtered triangle waveforms')

ax[1].plot(t, dx, 'b')
ax[1].plot(t, dy, 'g')
ax[1].set_ylabel('dx/dt')
ax[1].grid(True)

ax[1].set_xlabel('t')
plt.show()
