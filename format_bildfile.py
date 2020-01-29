#!/usr/bin/python

# make a bild file look better for making figures = change the size of spheres, make vectors into thicker cylinders


import sys

try:
    bildfile = open(sys.argv[1],'r').readlines()
    cyl_size = float(sys.argv[2])
    axis_size = float(sys.argv[3])
    sphere_size = float(sys.argv[4])
except:
    sys.exit('USAGE: bildfile_figure.py <bild file> <movement cylinder size> <axis cylinder size> <sphere size>\n use axis cylinder size of 0 to leave the axis markers as vectors')

output = open('FQ_{0}'.format(sys.argv[1].split('/')[-1]),'w')
for i in bildfile:
    line = i.split()
    if line[0] == '.cylinder':
        line = '.cylinder {0} {1}'.format(' '.join(line[1:-1]),cyl_size)
    elif line[0] == '.v' and axis_size !=0:
        line = '.cylinder {0} {1}'.format(' '.join(line[1:]),axis_size)
    elif line[0] == '.sphere':
        line = '{0} {1}'.format(' '.join(line[:-1]),sphere_size)
    else:
        line = ' '.join(line)
    output.write('{0}\n'.format(line))
output.close()
