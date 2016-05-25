#! /usr/bin/env python

import numpy
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import os
import sys

import seaborn

def ParseVariableBinaryHeader(header):

    header_detect = '#'

    header = header.strip("\n").split()

    assert(header[0] == header_detect)

    name = header[1]
    dtype = header[2]
    nb_components = int(header[3])
    nx = int(header[4])
    ny = int(header[5])
    nz = int(header[6])

    return name, dtype, nb_components, nx, ny, nz

receiver_filename = "output/receivers.dat"
receiver_file = open(receiver_filename, 'r')

readlines = receiver_file.read().split("\n")

receiver_file.close()

temp_filename = "tmp.binary"
tempfile = open(temp_filename, 'wb')

# Parse header
header = readlines[0]
name, dtype, nb_components, nx, ny, nz = ParseVariableBinaryHeader(header)

# Write data without header
for line in readlines[1:]:
    tempfile.write(line + "\n")

tempfile.close()
tempfile = open(temp_filename, 'rb')

data = numpy.fromfile(tempfile, dtype = 'float_')

tempfile.close()

if os.path.exists(temp_filename):
    print "Removing temporary file " + str(temp_filename)
    os.remove(temp_filename)

print data.shape, nx, nz

data = data.reshape(nz, nx)

amplitude_max = max(numpy.amax(data), - numpy.amin(data))

print "amplitude_max=", amplitude_max

rcv_ids = {'rcv1': (nx / 2, 'blue'),
           'rcv2': (nx / 4, 'red'),
           # 'rcv3': (3 * nx / 5, 'green'),
           # 'rcv4': (1200, 'orange'),
           # 'rcv5': (800, 'purple'),
       }

with seaborn.axes_style("dark"):

    fig, (ax1, ax2) = plt.subplots(1, 2)

    cmap = 'gray'
    ax1.imshow(data, cmap = cmap, interpolation = 'none', aspect = 'auto', vmin = - 0.04 * amplitude_max, vmax = 0.04 * amplitude_max)

    for key, value in rcv_ids.iteritems():
        rcv_id, color = value
        ax1.plot([rcv_id, rcv_id], [0.0, nz], color = color, linewidth = 2)
    ax1.set_xlim([0,nx])
    ax1.set_ylim([nz,0])

    cnt = 1

    for key, value in rcv_ids.iteritems():
        rcv_id, color = value
        offset = numpy.power(-1.0, cnt) * (2.0 * amplitude_max) * (cnt / 2)
        ax2.plot(offset + data[:, rcv_id], color = color, linewidth = 2, label = key)
        cnt += 1
    plt.legend()

# fig.set_size_inches(8, 4)
# fig.savefig('./test.png', dpi=100)

plt.show()

    


