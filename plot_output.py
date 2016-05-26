#! /usr/bin/env python

import numpy
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import optparse
import os
import sys

import seaborn

HEADER_DETECT = '#'

def ParseVariableBinaryMasterHeader(header):

    header = header.strip("\n").split()
    assert(header[0] == HEADER_DETECT)
    
    nb_variables = int(header[1])

    return nb_variables
    

def ParseVariableBinaryHeader(header):

    header = header.strip("\n").split()
    assert(header[0] == HEADER_DETECT)

    name = header[1]
    dtype = header[2]
    nb_components = int(header[3])
    nx = int(header[4])
    ny = int(header[5])
    nz = int(header[6])
    xmin = float(header[7])
    xmax = float(header[8])
    ymin = float(header[9])
    ymax = float(header[10])
    zmin = float(header[11])
    zmax = float(header[12])

    return name, dtype, nb_components, nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax

parser = optparse.OptionParser(usage="usage: %prog filename")
(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("wrong number of arguments")

out_filename = args[0]
out_file = open(out_filename, 'rb')

lines = out_file.read().split("\n")

nb_variables = ParseVariableBinaryMasterHeader(lines[0])
print "Found " + str(nb_variables) + " variables in file"
variable_names = []

for i in range(nb_variables):
    variable_header = lines[1 + i]
    print variable_header.strip("\n")
    name, dtype, nb_components, nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax = ParseVariableBinaryHeader(variable_header)
    variable_names.append(name)

print variable_names
    
out_file.close()

temp_filename = "tmp.binary"

if os.path.exists(temp_filename):
    os.remove(temp_filename)

tempfile = open(temp_filename, 'wb')
for line in lines[1 + nb_variables:]:
    tempfile.write(line + "\n")
tempfile.close()

tempfile = open(temp_filename, 'rb')

data = numpy.fromfile(tempfile, dtype = 'float_')
data = data.reshape(nb_variables, nz, nx)

tempfile.close()

variables = {}

for i in range(nb_variables):
    print "Adding variable " + str(variable_names[i]) + " to data"
    variables[variable_names[i]] = data[i,:,:]

variable_name_to_plot = "pressure_0"

amplitude_max = max(numpy.amax(variables[variable_name_to_plot]), - numpy.amin(variables[variable_name_to_plot]))

print "amplitude_max=", amplitude_max

line_ids = {'line1': (nx / 2, 'blue'),
           'line2': (nx / 4, 'orange')}

with seaborn.axes_style("dark"):

    fig, (ax1, ax2) = plt.subplots(1, 2)

    cmap = 'gray'
    aspect_ratio = ((zmax - zmin) * nx) / ((xmax - xmin) * nz)
    ax1.imshow(variables[variable_name_to_plot], cmap = cmap, vmin = - 0.1 * amplitude_max, vmax = 0.1 * amplitude_max, aspect = aspect_ratio)

    for key, value in line_ids.iteritems():
        line_id, color = value
        ax1.plot([line_id, line_id], [0.0, nz], color = color, linewidth = 2)
        ax1.set_xlim([0,nx])
        ax1.set_ylim([nz,0])
        ax1.set_title("Snapshot of pressure at final time")
    variable_names_1D_plot = ["pressure_0"]

    cnt = 1

    for key, value in line_ids.iteritems():
        line_id, color = value
        offset = numpy.power(-1.0, cnt) * (2.0 * amplitude_max) * (cnt / 2)
        for var in variable_names_1D_plot:
            if var == "pressure_ref":
                color = 'red'
            ax2.plot(offset + variables[var][:, line_id], color = color, linewidth = 1, label = var)
            ax2.scatter(range(nz), offset + variables[var][:, line_id], color = color, linewidth = 1)

        cnt += 1
    ax2.set_title("Snapshot of pressure at final time on 2 extracted lines")
    plt.legend()


# fig.set_size_inches(8, 4)
# fig.savefig('./test.png', dpi=100)

plt.show()

    


