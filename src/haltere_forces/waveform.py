import numpy as np
import scipy.signal as signal


def triangle(t, amplitude=1.0, period=1.0):
    """
    Generates a triangle waveform evaluated at the specified time points. 

    Parameters
    ----------
    t:  array like
        Input array of time points
    amplitude: float
        The amplitude (peak value) of the triangle waveform
    period: float
        The period of triangle waveform

    Returns
    ------
    x: array like
       Output array of triangle waveform values evaluated at times t.
    
    """
    s = ((t + 0.25*period) % period)/period
    x = amplitude*(np.where(s > 0.5, 1.0-s, s) - 0.25)/0.25
    return x 


def filtered_triangle(num_pt, num_cycle=1, amplitude=1.0, period=1.0, 
        cutoff_frequency=1.0, rescale=True):
    """
    Generates a lowpass filtered (zero phase delay)  triangle waveform. 

    Parameters
    ----------
    num_pt: int
        The number of points in the waveform
    num_cycle: int 
        The number of cycles in waveform
    amplitude: float
        The amplitude (peak value) of the triangle waveform
    period: float
        The period of triangle waveform
    cutoff_frequency: float
        the cutoff frequency of the lowpass filter
    rescale: bool
        whether or not to rescale the amplitude of the waveform
        so that it is equal to amplitude after filtering (default=True). 

    Returns
    ------
    t: array like
       Output array of time values for triangle waveform

    x: array like
       Output array of triangle waveform values evaluated at times t.

    """

    if not type(num_cycle) == int and num_cycle > 0:
        raise ValueError('num_cycle must be an integer > 0')

    # Create time points with padding (pre and nxt) to remove end effects. 
    t_mid = np.linspace(0,num_cycle*period, num_pt)
    dt = t_mid[1] - t_mid[0]
    t_pre = t_mid - t_mid[-1] - dt
    t_nxt = t_mid + t_mid[-1] + dt
    t_all = np.hstack((t_pre, t_mid, t_nxt))

    # Create Triangle waveform and filter it with forward backward filter. 
    x_all = triangle(t_all, amplitude, period)
    b, a  = signal.butter(1, 2*cutoff_frequency, btype='lowpass', output='ba', fs=1/dt)
    x_filt_all = signal.filtfilt(b, a, x_all)
    x_filt_mid = x_filt_all[num_pt:2*num_pt]

    if rescale:
        x_filt_max = np.absolute(x_filt_mid).max()
        x_filt_mid = amplitude*x_filt_mid/x_filt_max
    return t_mid, x_filt_mid



    





