# Copyright (C) 2015-2016: The University of Edinburgh
#                 Authors: Craig Warren and Antonis Giannopoulos
#
# This file is part of gprMax.
#
# gprMax is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# gprMax is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with gprMax.  If not, see <http://www.gnu.org/licenses/>.

import h5py

from gprMax._version import __version__


def write_hdf5(outputfile, G):
    """Write an output file in HDF5 format.

    Args:
        outputfile (str): Name of the output file.
        G (class): Grid class instance - holds essential parameters describing the model.
    """

    f = h5py.File(outputfile, 'w')
    f.attrs['gprMax'] = __version__
    f.attrs['Title'] = G.title
    f.attrs['Iterations'] = G.iterations
    f.attrs['nx, ny, nz'] = (G.nx, G.ny, G.nz)
    f.attrs['dx, dy, dz'] = (G.dx, G.dy, G.dz)
    f.attrs['dt'] = G.dt
    nsrc = len(G.voltagesources + G.hertziandipoles + G.magneticdipoles + G.transmissionlines)
    f.attrs['nsrc'] = nsrc
    f.attrs['nrx'] = len(G.rxs)
    f.attrs['srcsteps'] = (G.srcstepx, G.srcstepy, G.srcstepz)
    f.attrs['rxsteps'] = (G.rxstepx, G.rxstepy, G.rxstepz)

    # Create group for sources (except transmission lines); add type and positional data attributes
    srclist = G.voltagesources + G.hertziandipoles + G.magneticdipoles
    for srcindex, src in enumerate(srclist):
        grp = f.create_group('/srcs/src{}'.format(srcindex + 1))
        grp.attrs['Type'] = type(src).__name__
        grp.attrs['Position'] = (src.xcoord * G.dx, src.ycoord * G.dy, src.zcoord * G.dz)

    # Create group for transmission lines; add positional data, line resistance and line discretisation attributes; write arrays for line voltages and currents
    if G.transmissionlines:
        for tlindex, tl in enumerate(G.transmissionlines):
            grp = f.create_group('/tls/tl{}'.format(tlindex + 1))
            grp.attrs['Position'] = (tl.xcoord * G.dx, tl.ycoord * G.dy, tl.zcoord * G.dz)
            grp.attrs['Resistance'] = tl.resistance
            grp.attrs['dl'] = tl.dl
            # Save incident voltage and current
            grp['Vinc'] = tl.Vinc
            grp['Iinc'] = tl.Iinc
            # Save total voltage and current
            grp['Vtotal'] = tl.Vtotal
            grp['Itotal'] = tl.Itotal

            # Save the total voltage and current if available
            try:
                grp['V'] = tl.V
            except AttributeError:
                pass

            try:
                grp['I'] = tl.I
            except AttributeError:
                pass

    # Create group, add positional data and write field component arrays for receivers
    for rxindex, rx in enumerate(G.rxs):
        grp = f.create_group('/rxs/rx{}'.format(rxindex + 1))
        if rx.ID:
            grp.attrs['Name'] = rx.ID
        grp.attrs['Position'] = (rx.xcoord * G.dx, rx.ycoord * G.dy, rx.zcoord * G.dz)

        # Copy field values into hdf5 file
        fields = ['Ex', 'Ey', 'Ez', 'Hx', 'Hy', 'Hz', 'Ix', 'Iy', 'Iz']

        for field in fields:
            if field in rx.outputs:
                grp[field] = rx.outputs[field]
