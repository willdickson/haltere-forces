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
        try:
            shift = self.param['shift']
        except KeyError:
            shift = 0.25

        match self.param['waveform']:
            case 'triangle':
                angles = waveform.triangle(self.t, amplitude, period, shift=shift) 
            case 'filtered_triangle':
                num_pt = self.param['num_pt']
                num_cycle = self.param['num_cycle']
                cutoff_freq = self.param['cutoff_freq']
                _, angles = waveform.filtered_triangle(
                        num_pt=num_pt, 
                        num_cycle=num_cycle, 
                        amplitude=amplitude, 
                        period=period,
                        shift=shift, 
                        cutoff_frequency=cutoff_freq,
                        )
            case _:
                raise ValueError(f'unknown waveform, {self.param["waveform"]}')
        return angles

    @property
    def axis_left(self):
        num_pt = self.param['num_pt']
        axis = np.array([-1.0, 0.0, 0.0])
        rot_axis_angle = np.zeros((num_pt, 3))
        rot_axis_angle[:,1] = -self.angle
        qrot_flap = qn.array.from_axis_angle(rot_axis_angle)
        axis = qrot_flap.rotate(axis)
        return axis

    @property
    def axis_right(self):
        num_pt = self.param['num_pt']
        axis = np.array([ 1.0, 0.0, 0.0])
        rot_axis_angle = np.zeros((num_pt, 3))
        rot_axis_angle[:,1] = self.angle
        qrot_flap = qn.array.from_axis_angle(rot_axis_angle)
        axis = qrot_flap.rotate(axis)
        return axis

    @property
    def rot_axis_left(self): 
        tilt_angle = self.param['tilt_angle']
        axis = np.array([0.0, 1.0, 0.0])
        qrot_tilt = qn.array.from_axis_angle([0.0, 0.0, tilt_angle])
        axis = qrot_tilt.rotate(axis)
        return axis
        
    @property 
    def rot_axis_right(self): 
        tilt_angle = self.param['tilt_angle']
        axis = np.array([0.0, 1.0, 0.0])
        qrot_tilt = qn.array.from_axis_angle([0.0, 0.0, -tilt_angle])
        axis = qrot_tilt.rotate(axis)
        return axis

    @property 
    def pos_left(self):
        length = self.param['length']
        separation = self.param['separation']
        tilt_angle = self.param['tilt_angle']
        pos = self.axis_left*length
        qrot_tilt = qn.array.from_axis_angle([0.0, 0.0, tilt_angle])
        pos = qrot_tilt.rotate(pos)
        pos -= np.array([0.5*separation, 0.0, 0.0])
        return pos

    @property
    def pos_right(self):
        length = self.param['length']
        separation = self.param['separation']
        tilt_angle = self.param['tilt_angle']
        pos = self.axis_right*length
        qrot_tilt = qn.array.from_axis_angle([0.0, 0.0, -tilt_angle])
        pos = qrot_tilt.rotate(pos)
        pos += np.array([0.5*separation, 0.0, 0.0])
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

    def force_left(self, omega, domega=None, lin_acc=None):
        force = calc_haltere_force(
                self.param['mass'], 
                self.pos_left, 
                self.vel_left,
                self.acc_left, 
                self.axis_left, 
                self.rot_axis_left, 
                omega, 
                domega,
                lin_acc
                )
        return force

    def force_right(self, omega, domega=None, lin_acc=None):
        force = calc_haltere_force(
                self.param['mass'], 
                self.pos_right, 
                self.vel_right,
                self.acc_right, 
                self.axis_right,
                self.rot_axis_right, 
                omega, 
                domega,
                lin_acc
                )
        return force

    def force(self, omega, domega=None, lin_acc=None):
        force = {
                'left'  : self.force_left(omega, domega, lin_acc), 
                'right' : self.force_right(omega, domega, lin_acc), 
                }
        return force

    
# -----------------------------------------------------------------------------

def calc_haltere_force(m, h_pos, h_vel, h_acc, h_axis, h_rot_axis, omega, 
        domega=None, lin_acc=None):
    """
    Computes the forces on a haltere given the mass, the position, velocity
    and acceleration vectors of the haltere, the angular velocity and angular
    acceleration of the body and the body linear acceleration.

    Parameters:
    m : float
        the mass of the haltere
    h_pos : array_like
        shape (N,3) array of haltere position vectors
    h_vel : array_like
        shape (N,3) array of haltere velocity vectors
    h_acc : array_like
        shape (N,3) array of haltere acceleration vectors
    h_axis : array_like
        haltere stalk axis (3,1) array
    h_rot_axis : array_like
        haltere rotation axis (3,) array
    omega : array_like
        shape (3,), (1,3) or (N,3)  array of body angular velocities
    domega : array_like
        shape (3,), (1,3) or (N,3) array of body angular accelerations 
    lin_acc: array_like
        shape (3,), (1,3) or (N,3) array of body linear accelerations

    Returns:
    force : dict
        dictionary of (N,3) force arrays
        force = {
            'gravity'  : gravitational force, 
            'primary'  : force due to haltere acceleration,
            'lin_acc'  : force due to fly's linear acceleration,
            'ang_acc'  : force due to fly's angular acceleration,
            'centrif'  : centrifugal force on haltere,
            'coriolis  : coriolis forces on haltere,
        }

    """
    # Gravitational acceleration vector
    g = m*np.array([0.0, 0.0, -sp.constants.g])

    # Get array size and check shapes of h_pos, h_vel and h_acc
    n = h_pos.shape[0]
    check_shape(h_pos, (n, 3))
    check_shape(h_vel, (n, 3))
    check_shape(h_acc, (n, 3))

    # If domega or lin_acc aren't given assume they are zeros
    if domega is None: 
        domega = np.zeros_like(omega) 
    if lin_acc is None: 
        lin_acc = np.zeros_like(omega)

    # Reshape omega, domega, and lin_acc to (n,3) in necessary
    _omega  = reshape_to_nx3(n, omega) 
    _domega = reshape_to_nx3(n, domega)
    _lin_acc = reshape_to_nx3(n, lin_acc)
    _h_axis = reshape_to_nx3(n, h_axis)
    _h_rot_axis = reshape_to_nx3(n, h_rot_axis)

    ## Compute the forces
    f_gravity  =  m*g
    f_primary  = -m*h_acc 
    f_lin_acc  = -m*lin_acc
    f_ang_acc  = -m*np.cross(_domega, h_pos)
    f_centrif  = -m*np.cross(_omega, np.cross(_omega, h_pos))
    f_coriolis = -2.0*m*np.cross(_omega, h_vel)
    f_total = f_gravity + f_primary + f_lin_acc + f_ang_acc + f_centrif + f_coriolis 

    #print(_omega[:5,:])
    #print(h_vel[:5,:])
    #print(f_coriolis[:5,:])
    #print()

    force = {
            'total'       : f_total,
            'gravity'     : f_gravity,  
            'primary'     : f_primary, 
            'linear_acc'  : f_lin_acc, 
            'angular_acc' : f_ang_acc, 
            'centrifugal' : f_centrif, 
            'coriolis'    : f_coriolis, 
            }
    force_comp_names = list(force.keys())
    force['radial']  = {}
    force['lateral'] = {}
    for k in force_comp_names:
        force['radial'][k]  = project(force[k], _h_axis)
        force['lateral'][k] = project(force[k], _h_rot_axis)
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



def project(a,b):
    """
    Get projection of vector array a onto vectory array b

    Parameters
    a : array_like
        (n,k) array of vectors where each a[i,:] is length k vector
    b : array_like
        (n,k) array of vectory where each b[i,:] is length k vector

    Returns
    p : array_like
        (n,) array of projections of a vectors onto b vectors
    """
    b_norm = np.linalg.norm(b,axis=1)
    b_unit = b/b_norm[:, np.newaxis]
    p = np.sum(a*b_unit, axis=1)
    return p

   

    
















