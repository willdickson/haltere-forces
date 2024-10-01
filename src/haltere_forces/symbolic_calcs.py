import sympy

sympy.init_printing()


def simple_pitch_calc():

    t = sympy.Symbol('t')                  # time
    L = sympy.Symbol('L')                  # haltere length
    b = sympy.Symbol('b')                  # haltere base separation
    m = sympy.Symbol('m')                  # haltere mass
    beta = sympy.Symbol('beta')            # haltere tilt angle
    omega_x = sympy.Symbol('omega_x')      # fly x-component of angular velocity
    theta = sympy.Function('theta')(t)     # haltere angular position function

    # Fly's angular velocity vector
    omega = sympy.Matrix([omega_x, 0, 0])

    # Right and left haltere rotation matries for caudad tilt
    r_rot_mat = sympy.rot_ccw_axis3(-beta)
    l_rot_mat = sympy.rot_ccw_axis3( beta)

    # Right haltere position, velocity and rotation axis vectors
    r_pos = sympy.Matrix([L*sympy.cos(theta), 0, -L*sympy.sin(theta)])
    r_pos = r_rot_mat*r_pos
    r_pos = r_pos + sympy.Matrix([b/2, 0, 0])
    r_vel = sympy.diff(r_pos, t)
    r_axs = sympy.Matrix([0, 1, 0])
    r_axs = r_rot_mat*r_axs

    # Left haltere position and velocity vectors
    l_pos = sympy.Matrix([-L*sympy.cos(theta), 0, -L*sympy.sin(theta)])
    l_pos = l_rot_mat*l_pos
    l_pos = l_pos - sympy.Matrix([b/2, 0, 0])
    l_vel = sympy.diff(l_pos, t)
    l_axs = sympy.Matrix([0, 1, 0])
    l_axs = l_rot_mat*l_axs

    # Coriolis forces for right and left wings 
    r_coriolis = 2*m*omega.cross(r_vel)
    l_coriolis = 2*m*omega.cross(l_vel)

    r_coriolis_lateral = r_coriolis.dot(r_axs)
    l_coriolis_lateral = l_coriolis.dot(l_axs)

    print()
    print('haltere angle = ')
    sympy.pprint(theta)

    print()
    print('right haltere tilt rotation matrix')
    sympy.pprint(r_rot_mat)

    print()
    print('left haltere tilt rotation matrix')
    sympy.pprint(l_rot_mat)

    print()
    print('fly angular velocity = ')
    sympy.pprint(omega)

    print()
    print('right haltere rotation axis')
    sympy.pprint(r_axs)

    print()
    print('left haltere rotation axis')
    sympy.pprint(l_axs)
    
    print()
    print('right haltere position =')
    sympy.pprint(r_pos)

    print()
    print('left haltere position =')
    sympy.pprint(l_pos)

    print()
    print('right haltere velocity =')
    sympy.pprint(r_vel)

    print()
    print('left haltere velocity =')
    sympy.pprint(l_vel)

    print()
    print('right coriolis force = ')
    sympy.pprint(r_coriolis)

    print()
    print('left coriolis force = ')
    sympy.pprint(l_coriolis)

    print()
    print('right coriolis lateral = ')
    sympy.pprint(r_coriolis_lateral)

    print()
    print('left coriolis lateral = ')
    sympy.pprint(l_coriolis_lateral)







