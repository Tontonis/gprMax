from gprMax import GPRObjectCreator

oc = GPRObjectCreator()


def create_dipole(x, y, z, dl=0.001):

    centre = [x, y, z]

    # dipole length
    length = 0.095

    # Source excitation and type
    f = 3e9
    waveform = oc.create('waveform', 'gaussian', 1, f, 'mypulse')
    transmission_line = oc.create('transmission_line', 'z', centre[0], centre[1], centre[2] - dl, 50, 'mypulse')

    # dipole edge
    x1 = centre[0]
    y1 = centre[1]
    z1 = centre[2] - length / 2
    z1 = centre[2] - length / 2

    x2 = centre[0]
    y2 = centre[1]
    z2 = centre[2] + length / 2

    pec_edge = oc.create('edge', x1, y1, z1, x2, y2, z2, 'pec')

    # dipole gap
    x3 = centre[0]
    y3 = centre[1]
    z3 = centre[2] - dl

    x4 = centre[0]
    y4 = centre[1]
    z4 = centre[2]

    free_edge = oc.create('edge', x3, y3, z3, x4, y4, z4, 'free_space')

    antenna = oc.create('wrapper')
    antenna.add(waveform, transmission_line, pec_edge, free_edge)
    return antenna
