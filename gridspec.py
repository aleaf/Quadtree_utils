__author__ = 'aleaf'

import pandas as pd
import numpy as np
from shapely.geometry import box

tsf = '/Users/aleaf/Documents/ATLData/Wbasin/grid02qtg.tsf'
nod_file = '/Users/aleaf/Documents/ATLData/Wbasin/output_usgdata/grid02qtg.nod'

class gridSPC:

    def __init__(self, nod_file):
        self.nod_file = nod_file
        self.outfile = nod_file[:-4] + '.spc'

    def read_nod(self):

        self.nods = pd.read_csv(self.nod_file, skiprows=1,
                           delim_whitespace=True,
                           names=['node', 'layer', 'centroid_x', 'centroid_y', 'centroid_z', 'delta_x', 'delta_y', 'delta_z'],
                           header=None)

        self.nnode = len(self.nods)
        self.nlay = len(np.unique(self.nods['layer']))
        self.iz = 1
        self.ic = 1


        self.nvertex = self.nnode * 8 # for now, write all eight vertices for each node

        # build a list of vertices for each node

        def vertices(x):
            xmin, xmax = x['centroid_x'] - 0.5 * x['delta_x'], x['centroid_x'] + 0.5 * x['delta_x']
            ymin, ymax = x['centroid_y'] - 0.5 * x['delta_y'], x['centroid_y'] + 0.5 * x['delta_y']
            zmin, zmax = x['centroid_z'] - 0.5 * x['delta_z'], x['centroid_z'] + 0.5 * x['delta_z']

            x['verts'] = [(xmin, ymin, zmin),
                          (xmin, ymax, zmin),
                          (xmax, ymin, zmin),
                          (xmax, ymax, zmin),
                          (xmin, ymin, zmax),
                          (xmin, ymax, zmax),
                          (xmax, ymin, zmax),
                          (xmax, ymax, zmax)]


        self.nods.apply(vertices, axis=1)
        #verts = [vertices(r) for i, r in nods.iterrows()]


    def write_byline(self, nlay, iz=1, ic=1):

        nods = open(self.nod_file, 'r')

        nnods = int(nods.readline().strip().split()[0])
        nvertex = nnods * 8

        ofp = open(self.outfile, 'w')
        ofp.write('UNSTRUCTURED GWF\n')
        ofp.write('{} {} {} {}\n'.format(nnods, nlay, iz, ic))
        ofp.write('{}\n'.format(nvertex))


        def vertices(line):

            cx, cy, cz, dx, dy, dz = map(float, line.split()[2:8])

            xmin, xmax = cx - 0.5 * dx, cx + 0.5 * dx
            ymin, ymax = cy - 0.5 * dy, cy + 0.5 * dy
            zmin, zmax = cz - 0.5 * dz, cz + 0.5 * dz

            verts = [(xmin, ymin, zmin),
                     (xmin, ymax, zmin),
                     (xmax, ymin, zmin),
                     (xmax, ymax, zmin),
                     (xmin, ymin, zmax),
                     (xmin, ymax, zmax),
                     (xmax, ymin, zmax),
                     (xmax, ymax, zmax)]

            return verts


        for line in nods:

            verts = vertices(line)
            [[ofp.write(' '.join(map(str, v)) + '\n')] for v in verts]

        i = 1
        inode = 1
        for line in nods:

            lay = int(line.split()[1])
            x, y, z = map(float, line.split()[2:5])

            ofp.write('{:.0f} {:.0f} {:.0f} {:.0f} {:.0f} 8'.format(inode, x, y, z, lay))

            ofp.write([' {:.0f}'.format(c) for c in np.arange(i, i + 8)])
            ofp.write('\n')
            i += 8
            inode += 1



    def write(self):

        ofp = open(self.outfile, 'w')

        ofp.write('UNSTRUCTURED GWF\n')
        ofp.write('{} {} {} {}\n'.format(self.nnode, self.nlay, self.iz, self.ic))
        ofp.write('{}\n'.format(self.nvertex))

        # write out list of vertices (including duplicates)
        [ofp.write(' '.join(v) + '\n') for v in self.nods['verts']]

        # write out node information
        i = 1
        for i, r in self.nods.iterrows():

            ofp.write('{:.0f} {:.0f} {:.0f} {:.0f} {:.0f} {:.0f}'.format())

            ofp.write([' {:.0f}'.format(c) for c in np.arange(i, i + 8)])
            ofp.write('\n')
            i += 8

        ofp.close()

