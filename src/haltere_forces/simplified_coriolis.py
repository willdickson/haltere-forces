import numpy as np

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
    num_pts = len(theta)
    return 2.0*length*mass*omega_x*dtheta_dt*np.cos(theta)*np.cos(beta)



