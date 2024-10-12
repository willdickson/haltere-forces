import numpy as np


def coriolis_from_yaw(mass, length, beta, omega_z, theta, dtheta_dt, side):
    """
    Calculates the coriolis force due to a constant yaw rotation (rotation 
    about the z-axis). 

    Parameters
    ----------
    mass : float
        the mass of the haltere end knob 
    length : float
        the length of the haltere, from point of rotation to center of end knob.
    beta : float
        haltere tilt angle
    omega_z : float
        pitch rate about the z-axis 
    theta : array_like
        haltere position angles (rad). 
    dtheta_dt : array_like
        haltere angular velocities (rad/sec)
    side: string
        haltere side, 'left' or 'right'

    Returns
    -------
    f : array_like
        coriolis forces 

    """
    num_pts = len(theta)
    f = np.zeros((num_pts, 3))
    sing_y = 1 if side == 'right' else -1
    f[:,0] = 2.0*length*mass*omega_z*dtheta_dt*np.sin(theta)*np.sin(beta)
    f[:,1] = sign_y*2.0*length*mass*omega_z*dtheta_dt*np.sin(theta)*np.cos(beta)
    return f


def lateral_coriolis_from_yaw(mass, length, omega_z, theta, dtheta_dt, side):
    """
    Calculates the lateral component of the coriolis force due to a constant
    yaw rotation (rotation about the z-axis). 

    Parameters
    ----------
    mass : float
        the mass of the haltere end knob 
    length : float
        the length of the haltere, from point of rotation to center of end knob.
    omega_z : float
        pitch rate about the z-axis 
    theta : array_like
        haltere position angles (rad). 
    dtheta_dt : array_like
        haltere angular velocities (rad/sec)
    side: string
        haltere side, 'left' or 'right'

    Returns
    -------
    f : array_like
        coriolis forces 

    """
    sign = 1 if side == 'right' else -1
    return 2.0*sign*length*mass*omega_z*dtheta_dt*np.sin(theta)


def coriolis_from_pitch(mass, length, beta, omega_x, theta, dtheta_dt):
    """
    Calculates the coriolis force due to a constant pitch rotation (rotation
    about the x-axis).

    Parameters
    ----------
    mass : float
        the mass of the haltere end knob 
    length : float
        the length of the haltere, from point of rotation to center of end knob.
    beta : float
        haltere tilt angle
    omega_x : float
        pitch rate about the x-axis 
    theta : array_like
        haltere position angles (rad). 
    dtheta_dt : array_like
        haltere angular velocities (rad/sec)

    Returns
    -------
    f : array_like
        coriolis forces 
        
    """
    num_pts = len(theta)
    f = np.zeros((num_pts, 3))
    f[:,1] = 2.0*length*mass*omega_x*dtheta_dt*np.cos(theta)
    f[:,2] = -2.0*length*mass*omega_x*dtheta_dt*np.sin(theta)*np.sin(beta)
    return f


def lateral_coriolis_from_pitch(mass, length, beta, omega_x, theta, dtheta_dt):
    """
    Calculates the lateral component of the coriolis force due to a constant
    pitch rotation (rotation about the x-axis). 

    Parameters
    ----------
    mass : float
        the mass of the haltere end knob 
    length : float
        the length of the haltere, from point of rotation to center of end knob.
    beta : float
        haltere tilt angle
    omega_x : float
        pitch rate about the x-axis 
    theta : array_like
        haltere position angles (rad). 
    dtheta_dt : array_like
        haltere angular velocities (rad/sec)

    Returns
    -------
    f : array_like
        lateral component of coriolis forces   
        
    """
    return 2.0*length*mass*omega_x*dtheta_dt*np.cos(theta)*np.cos(beta)


def coriolis_from_roll(mass, length, beta, omega_y, theta, dtheta_dt, side):
    """
    Calculates the coriolis force due to a constant roll rotation (rotation 
    about the y-axis). 

    Parameters
    ----------
    mass : float
        the mass of the haltere end knob 
    length : float
        the length of the haltere, from point of rotation to center of end knob.
    beta : float
        haltere tilt angle
    omega_y : float
        pitch rate about the y-axis 
    theta : array_like
        haltere position angles (rad). 
    dtheta_dt : array_like
        haltere angular velocities (rad/sec)
    side: string
        haltere side, 'left' or 'right'

    Returns
    -------
    f : array_like
        coriolis forces 

    """
    num_pts = len(theta)
    f = np.zeros((num_pts, 3))
    sign_z = -1 if side == 'right' else 1
    f[:,0] = -2.0*length*mass*omega_y*dtheta_dt*np.cos(theta)*np.sin(beta)
    f[:,2] = sign_z*2.0*length*mass*omega_y*dtheta_dt*np.cos(theta)*np.sin(beta)
    return f


def lateral_coriolis_from_roll(mass, length, beta, omega_y, theta, dtheta_dt, side):
    """
    Calculates lateral component of the coriolis force due to a constant roll
    rotation (rotation about the y-axis). 

    Parameters
    ----------
    mass : float
        the mass of the haltere end knob 
    length : float
        the length of the haltere, from point of rotation to center of end knob.
    beta : float
        haltere tilt angle
    omega_y : float
        pitch rate about the y-axis 
    theta : array_like
        haltere position angles (rad). 
    dtheta_dt : array_like
        haltere angular velocities (rad/sec)
    side: string
        haltere side, 'left' or 'right'

    Returns
    -------
    f : array_like
        coriolis forces 

    """

    sign = -1 if side == 'right' else 1
    return 2.0*sign*length*mass*omega_y*dtheta_dt*np.cos(theta)*np.sin(beta)


