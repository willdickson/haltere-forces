import math 

def dprint(d):
    for k, v in d.items():
        print(f'{k}: {v:1.3e}')

def ellipsoid_volume(d1, d2, d3):
    return math.pi*d1*d2*d3/6.0



# -------------------------------------------------------------------------

scaling = 1/13450 # from (igor units) to (m)
density = 1000.0  # kg/m**2

# Igors haltere data
igor_values = {
        'length'         : 4.57, 
        'separation'     : 7.65,
        'major_diameter' : 3.0,
        'minor_diameter' : 1.7,
        }
scaled_values = {k: v*scaling for k, v in igor_values.items()}

# Get mass of haltere end bulb
d1 = scaled_values['major_diameter']
d2 = scaled_values['minor_diameter']
d3 = scaled_values['minor_diameter']
bulb_volume = ellipsoid_volume(d1, d2, d3)
scaled_values['mass'] = density*bulb_volume

print()
print('igors values')
print('------------')
dprint(igor_values)
print()
print('scaled values')
print('-------------')
dprint(scaled_values)
print()





