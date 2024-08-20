import numpy as np
import scipy as sp
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
    def dt(self):
        dt = 1.0/self.param['frequency']*self.param['num_cycle']/(self.param['num_pt']-1)
        return dt

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
    def pos_left(self):
        num_pt = self.param['num_pt']
        tilt_angle = self.param['tilt_angle']
        pos = np.array([-1.0, 0.0, 0.0])
        # Rotate left haltere through flapping motion
        axis_angle = np.zeros((num_pt, 3))
        axis_angle[:,1] = -self.angle
        qrot_flap = qn.array.from_axis_angle(axis_angle)
        pos = qrot_flap.rotate(pos)
        # Rotate left haltere for tilt angle
        qrot_tilt = qn.array.from_axis_angle([0.0, 0.0, tilt_angle])
        pos = qrot_tilt.rotate(pos)
        return pos

    @property
    def pos_right(self):
        num_pt = self.param['num_pt']
        tilt_angle = self.param['tilt_angle']
        pos = np.array([ 1.0, 0.0, 0.0])
        # Rotate right haltere through flapping motion
        axis_angle = np.zeros((num_pt, 3))
        axis_angle[:,1] =  self.angle
        qrot_flap = qn.array.from_axis_angle(axis_angle)
        pos = qrot_flap.rotate(pos)
        # Rotate right haltere for tilt angle
        qrot_tilt = qn.array.from_axis_angle([0.0, 0.0, -tilt_angle])
        pos = qrot_tilt.rotate(pos)
        return pos

    @property 
    def vel_left(self):
        return np.gradient(self.pos_left, axis=0)/self.dt

    @property
    def vel_right(self):
        return np.gradient(self.pos_right, axis=0)/self.dt

    @property
    def acc_left(self):
        return np.gradient(self.vel_left, axis=0)/self.dt

    @property
    def acc_right(self):
        return np.gradient(self.vel_right, axis=0)/self.dt

    @property
    def kinematics_left(self):
        kinematics = { 
                'pos'   : self.pos_left, 
                'vel'   : self.vel_left, 
                'acc'   : self.acc_left, 
                }
        return kinematics

    @property
    def kinematics_right(self):
        kinematics = { 
                'pos'   : self.pos_right, 
                'vel'   : self.vel_right, 
                'acc'   : self.acc_right,
                }
        return kinematics

    @property
    def kinematics(self):
        kinematics = {
                'left'  : self.kinematics_left, 
                'right' : self.kinematics_right, 
                }
        return kinematics

    def force_left(self, omega, domega=None, linacc=None):
        force = calc_haltere_force(
                self.param['mass'], 
                self.pos_left, 
                self.vel_left,
                self.acc_left, 
                omega, 
                domega,
                linacc
                )
        return force

    def force_right(self, omega, domega=None, linacc=None):
        force = calc_haltere_force(
                self.param['mass'], 
                self.pos_right, 
                self.vel_right,
                self.acc_right, 
                omega, 
                domega,
                linacc
                )
        return force

    def force(self, omega, domega=None, linacc=None):
        force = {
                'left'  : self.force_left, 
                'right' : self.force_right, 
                }
        return force

    
# -------------------------------------------------------------------

def calc_haltere_force(m, hpos, hvel, hacc, omega, domega=None, linacc=None):
    """
    Computes the forces on a haltere given the mass, the position, velocity
    and acceleration vectors of the haltere, the angular velocity and angular
    acceleration of the body and the body linear acceleration.

    Parameters:
    m : float
        the mass of the haltere
    hpos : array_like
        shape (N,3) array of haltere position vectors
    hvel : array_like
        shape (N,3) array of haltere velocity vectors
    hacc : array_like
        shape (N,3) array of haltere acceleration vectors
    omega : array_like
        shape (3,), (1,3) or (N,3)  array of body angular velocities
    domega : array_like
        shape (3,), (1,3) or (N,3) array of body angular accelerations 
    linacc: array_like
        shape (3,), (1,3) or (N,3) array of body linear accelerations

    Returns:
    force : array_like
        shape (N,3) array of body forces

    """
    # Gravitational acceleration vector
    g = m*np.array([0.0, 0.0, -sp.constants.g])

    # Get array size and check shapes of hpos, hvel and hacc
    n = hpos.shape[0]
    check_shape(hpos, (n, 3))
    check_shape(hvel, (n, 3))
    check_shape(hacc, (n, 3))

    print(hpos)

    # If domega or linacc aren't given assume they are zeros
    if domega is None: 
        domega = np.zeros_like(omega) 
    if linacc is None: 
        linacc = np.zeros_like(omega)

    # Reshape omega, domega, and linacc to (n,3) in necessary
    _omega  = reshape_to_nx3(n, omega) 
    _domega = reshape_to_nx3(n, domega)
    _linacc = reshape_to_nx3(n, linacc)

    # Compute the forces
    force  = m*g - m*hacc - m*linacc 
    force -= m*np.cross(_domega, hpos)
    force -= m*np.cross(_omega, np.cross(_omega, hpos))
    force -= 2.0*m*np.cross(_omega, hvel)
    return force


def check_shape(a, shape):
    """
    Check to see if the shape of array is equal to shape. Raises
    ValueError is a.shape is not equal to shape.

    Parameters
    a : array_like
        the input array to be tested
    shape : the shape to test against

    """
    if a.shape != shape:
        raise ValueError(f'shape of array must be ({n},{m})')


def reshape_to_nx3(n, a):
    """
    Reshapes input array to shape (n,3) array. If input is (3,) or (1,3)
    the array is repeated n times. 

    Parameters 
    n : int
        number of rows in a
    a : array like
        the input array must be (3,), (1,3) or (n,3)

    Returns:
    _a : array_like
        the array a reshaped to be (n

    """
    match a.shape:
        case (3,):
            _a = np.repeat(a[np.newaxis, ...], n, axis=0)
        case (1,3): 
            _a = np.repeat(a, n, axis=0) 
        case (k,3) if k==n: 
            _a = a
        case _:
            raise ValueError('array, a, must be (3,), (1,3) or (n,3)')
    return _a



   

    
















