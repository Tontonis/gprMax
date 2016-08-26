import gprMax
from user_libs.dipole import create_dipole

scene = gprMax.Scene()
oc = gprMax.GPRObjectCreator()

dl = 0.001
pmlcells = 10
boundarycells = 40
domain_x = (2 * (pmlcells + boundarycells) * dl) + dl
domain_y = (2 * (pmlcells + boundarycells) * dl) + dl
domain_z = (2 * (pmlcells + boundarycells) * dl) + 0.1
centre = [domain_x / 2, domain_y / 2, domain_z / 2]

title = oc.create('title', 'dipole')
dxdydz = oc.create('dx_dy_dz', dl, dl, dl)
time_window = oc.create('time_window', 10e-9)
domain = oc.create('domain', domain_x, domain_y, domain_z)

antenna = create_dipole(centre[0], centre[1], centre[2])

scene.add(title, domain, dxdydz, time_window, antenna)

# Geometry view of entire domain
#geometry_view(0, 0, 0, domain[0], domain[1], domain[2], dxdydz[0], dxdydz[1], dxdydz[2], title, 'f')

gprMax.run(scene=scene)
