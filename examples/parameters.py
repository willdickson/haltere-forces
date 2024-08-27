import numpy as np

calliphora_param = {
        'num_pt'      : 100,    
        'num_cycle'   : 1, 
        'waveform'    : 'filtered_triangle',
        'cutoff_freq' : 2*130.0,             # (Hz) 
        'frequency'   : 130.0,               # (Hz)
        'amplitude'   : np.deg2rad(90.0),    # (deg)
        'mass'        : 5.89e-9,             # (kg) 
        'length'      : 1.2e-3,              # (m)
        'separation'  : 3.0e-3,              # (m)
        'tilt_angle'  : np.deg2rad(30.0),    # (rad)
        }

# Scaling parameters
length_calliphora = 17.0
length_drosophila = 3.0
length_sf = length_drosophila/length_calliphora
volume_sf = length_sf**3

drosophila_param = dict(calliphora_param)
drosophila_param['cutoff_freq'] = 2*200
drosophila_param['frequency']   = 200
drosophila_param['mass']        = calliphora_param['mass']*volume_sf 
drosophila_param['length']      = calliphora_param['length']*length_sf 
drosophila_param['separation']  = calliphora_param['separation']*length_sf 

