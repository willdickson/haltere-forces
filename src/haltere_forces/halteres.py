import numpy as np
import quaternionic as qn
from . import waveform

class Halteres:

    def __init__(self,param):
        self.param = param

    @property
    def t(self):
        period = 1.0/self.param['frequency']
        t = np.linspace(0, period*self.param['num_cycle'], self.param['num_pt']) 
        return t

    @property
    def angle(self):
        amplitude = self.param['amplitude']
        frequency = self.param['frequency']
        period = 1.0/frequency
        match self.param['waveform']:
            case 'triangle':
                angles = waveform.triangle(self.t, amplitude, period) 
            case 'filtered_triangle':
                num_pt = self.param['num_pt']
                num_cycle = self.param['num_cycle']
                cutoff_freq = self.param['cutoff_freq']
                _, angles = waveform.filtered_triangle(
                        num_pt=num_pt, 
                        num_cycle=num_cycle, 
                        amplitude=amplitude, 
                        period=period,
                        cutoff_frequency=cutoff_freq,
                        )
            case _:
                raise ValueError(f'unknown waveform, {self.param["waveform"]}')
        return angles


    @property
    def positions(self):
        num_pt = self.param['num_pt']
        tilt_angle = self.param['tilt_angle']

        # Left haltere positions
        pos_left = np.array([-1.0, 0.0, 0.0])

        # Rotate left haltere through flapping motion
        axis_angle_left = np.zeros((num_pt, 3))
        axis_angle_left[:,1] = -self.angle
        qrot_flap_left = qn.array.from_axis_angle(axis_angle_left)
        pos_left = qrot_flap_left.rotate(pos_left)

        # Rotate left haltere for tilt angle
        qrot_tilt_left = qn.array.from_axis_angle([0.0, 0.0, tilt_angle])
        pos_left = qrot_tilt_left.rotate(pos_left)

        # Right haltere positions
        pos_right = np.array([ 1.0, 0.0, 0.0])

        # Rotate right haltere through flapping motion
        axis_angle_right = np.zeros((num_pt, 3))
        axis_angle_right[:,1] =  self.angle
        qrot_flap_right = qn.array.from_axis_angle(axis_angle_right)
        pos_right = qrot_flap_right.rotate(pos_right)

        # Rotate right haltere for tilt angle
        qrot_tilt_right = qn.array.from_axis_angle([0.0, 0.0, -tilt_angle])
        pos_right = qrot_tilt_right.rotate(pos_right)

        return pos_left, pos_right













