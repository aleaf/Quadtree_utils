__author__ = 'aleaf'

import numpy as np


tsf = '/Users/aleaf/Documents/ATLData/Wbasin/grid02qtg.tsf'
nod_file = '/Users/aleaf/Documents/ATLData/Wbasin/output_usgdata/grid02qtg.nod'


class parentGrid(object):


    def __init__(self, nrow, ncol, nlay, xll=0, yll=0, dxy=1):

        self.nrow = nrow
        self.ncol = ncol
        self.nlay = nlay
        self.xll = xll
        self.yll = yll
        self.dxy = dxy


    def read_elevations(self, topfiles, botfiles):

        print "reading cell top and bottom elevations..."
        self.tops = np.array([], dtype=float)
        for i, f in enumerate(sorted(topfiles)):
            print f
            self.tops = np.append(self.tops, np.fromfile(f, sep=' '))

        self.bots = np.array([], dtype=float)
        for i, f in enumerate(sorted(botfiles)):
            print f
            self.bots = np.append(self.bots, np.fromfile(f, sep=' '))

        self.nnod = len(self.bots)



class qtCell:

    def __init__(self, r, c, qcode, nrow, dxy, xll=0, yll=0):

        # information on parent grid
        self.nrow = nrow
        self.dxy = dxy
        self.xll = xll
        self.yll = yll

        # information from tree structure file
        self.r = r
        self.c = c
        self.qcode = qcode

        # quadtree refinement information
        if self.qcode != 0:
            self.ref = len(self.qcode)
            self.nsgperbg = 2 ** self.ref
            self.sgdxy = self.dxy / self.nsgperbg

        # calculate vertices for parent cell
        self.x0 = self.dxy * self.c + self.xll
        self.x1 = self.dxy * (self.c - 1) + self.xll
        self.y0 = self.dxy * (self.nrow - self.r) + self.yll
        self.y1 = self.dxy * (self.nrow - (self.r - 1)) + self.yll


    def rangefinder(self):
        '''
        get row / column extent of refined cell within parent cell based on quadtree code
        '''
        isgmin = 0
        isgmax = isgmin + self.nsgperbg
        jsgmin = 0
        jsgmax = jsgmin + self.nsgperbg

        # get row column extent of refined cell based on code
        for q in self.qcode:
            if q is '1':  #upper left
                isgmax = isgmin + (isgmax - isgmin) / 2
                jsgmax = jsgmin + (jsgmax - jsgmin) / 2
            elif q is '2':  #upper right
                isgmax = isgmin + (isgmax - isgmin) / 2
                jsgmin = jsgmax - (jsgmax - jsgmin) / 2
            elif q is '3':  #lower left
                isgmin = isgmax - (isgmax - isgmin) / 2
                jsgmax = jsgmin + (jsgmax - jsgmin) / 2
            elif q is '4':  #lower right
                isgmin = isgmax - (isgmax - isgmin) / 2
                jsgmin = jsgmax - (jsgmax - jsgmin) / 2
            else:
                msg = 'Uknown quadrant code: ' + q + ' for cell ' + str(self.r, self.c)
                raise Exception(msg)

        return isgmin, isgmax, jsgmin, jsgmax


    def ts2coords(self):
        '''
        get vertices for refined cell using range finder and parent cell vertices
        '''

        # get the row / column extent of refine qt cell within parent cell
        isgmin, isgmax, jsgmin, jsgmax = self.rangefinder()

        # now convert refined cell row column extent to cell coordinates
        ymax, ymin, xmin, xmax = self.sgdxy * np.array([self.nsgperbg - isgmin,
                                                        self.nsgperbg - isgmax,
                                                        jsgmin,
                                                        jsgmax])
        # return cell coordinates as global coordinates
        return xmin + self.x0, ymin + self.y0, xmax + self.x0, ymax + self.y0


    def get_centroid(self, r, c):
        '''
        rcl (tuple) (row, column, layer)

        '''
        x = self.dxy * (c - 0.5)
        y = self.dxy * (self.nrow - r + 0.5)
        return x, y



class gridSPC(parentGrid):

    def write_from_nodfile(self, nod_file, iz=1, ic=1):
        '''
        write spc file from nod file (note: the node numbering is wrong compared to MODFLOW!
        '''
        nods = open(nod_file, 'r')

        nnods = int(nods.readline().strip().split()[0])
        nvertex = nnods * 8

        outfile = nod_file[:-4] + '.spc'
        ofp = open(outfile, 'w')
        ofp.write('UNSTRUCTURED GWF\n')
        ofp.write('{} {} {} {}\n'.format(nnods, self.nlay, iz, ic))
        ofp.write('{}\n'.format(nvertex))


        def vertices(line):

            cx, cy, cz, dx, dy, dz = map(float, line.split()[2:8])

            xmin, xmax = cx - 0.5 * dx, cx + 0.5 * dx
            ymin, ymax = cy - 0.5 * dy, cy + 0.5 * dy
            zmin, zmax = cz - 0.5 * dz, cz + 0.5 * dz

            verts = [(xmin, ymin, zmax),
                     (xmax, ymin, zmax),
                     (xmax, ymax, zmax),
                     (xmin, ymax, zmax),
                     (xmin, ymin, zmin),
                     (xmax, ymin, zmin),
                     (xmax, ymax, zmin),
                     (xmin, ymax, zmin)]

            return verts

        for line in nods:

            verts = vertices(line)
            [[ofp.write(' '.join(map(str, v)) + '\n')] for v in verts]

        nods = open(self.nod_file, 'r')
        nods.next()

        i = 1
        inode = 1
        for line in nods:

            lay = int(line.split()[1])
            x, y, z = map(float, line.split()[2:5])

            ofp.write('{:.0f} {:.2f} {:.2f} {:.2f} {:.0f} 8'.format(inode, x, y, z, lay))

            [ofp.write(' {:.0f}'.format(c)) for c in np.arange(i, i + 8)]
            ofp.write('\n')
            i += 8
            inode += 1

        ofp.close()

    def centroid(self, vertxy, vertz):
        '''
        return cell centroid from vertices
        '''
        xmin, ymin, = vertxy[0]
        zmax = vertz[0]
        xmax, ymax = vertxy[6]
        zmin = vertz[6]

        cx = 0.5 * (xmin + xmax)
        cy = 0.5 * (ymin + ymax)
        cz = 0.5 * (zmin + zmax)

        return cx, cy, cz


    def write_from_tsf(self, tsf_file, iz=1, ic=1, smooth=False, outfile=None):
        '''
        write grid spc file from tree structure file

        '''
        # total number of cell vertices to be listed
        self.nvertext = self.nnod * 8

        if not outfile:
            outfile = tsf_file[:-4] + '.spc'

        print 'reading quadtree grid information from {}...'.format(tsf_file)
        print 'writing output to {}...'.format(outfile)

        # open output file and write header
        ofp = open(outfile, 'w')
        ofp.write('UNSTRUCTURED GWF\n')
        ofp.write('{} {} {} {}\n'.format(self.nnod, self.nlay, iz, ic))
        ofp.write('{}\n'.format(self.nvertext))

        # open tsf file
        nods = open(tsf_file, 'r')
        nods.next() # skip the header

        allvertxy = np.zeros((self.nnod, 8), dtype=('float,float'))
        allvertz = np.zeros((self.nnod, 8), dtype=('float'))
        layer = np.zeros(self.nnod)
        for line in nods:

            inode = int(line.split(',')[0])

            # skip inactive nodes
            if inode == -1:
                continue

            r, c, l = map(int, line.split()[1].replace('(', '').replace(')', '').split(','))
            try:
                qcode = [q for q in line.strip().split()[2]]
            except:
                qcode = 0


            # initialize class for quadtree cell
            cell = qtCell(r, c, qcode, self.nrow, self.dxy, xll=self.xll, yll=self.yll)

            if qcode == 0:
                xmin, ymin, xmax, ymax = cell.x0, cell.y0, cell.x1, cell.y1

            else:
                # get coordinates of child cell relative to xmin/xmax of parent
                xmin, ymin, xmax, ymax = cell.ts2coords()

            # get zmin / zmax values
            zmin, zmax = self.tops[inode - 1], self.bots[inode - 1]

            # make list of cell vertices in ccw order, top then bottom
            verts = [(xmin, ymin, zmax),
                     (xmax, ymin, zmax),
                     (xmax, ymax, zmax),
                     (xmin, ymax, zmax),
                     (xmin, ymin, zmin),
                     (xmax, ymin, zmin),
                     (xmax, ymax, zmin),
                     (xmin, ymax, zmin)]

            vertsxy = [v[0:2] for v in verts]
            for i, v in enumerate(vertsxy):
                allvertxy[inode - 1, i] = v
                allvertz[inode - 1, i] = verts[i][2]
                layer[inode - 1] = l

        if smooth:

            #implement option to remove duplicate vertices here
            pass
        else:

            for n in range(len(allvertxy)):

                [ofp.write('{:.2f} {:.2f} {:.2f}\n'.format(allvertxy[n, v][0], allvertxy[n, v][1], allvertz[n, v]))
                 for v in range(8)]
                
            i = 1
            for n in range(len(allvertxy)):

                x, y, z = self.centroid(allvertxy[n], allvertz[n])

                ofp.write('{:.0f} {:.2f} {:.2f} {:.2f} {:.0f} 8'.format(n+1, x, y, z, layer[n]))

                [ofp.write(' {:.0f}'.format(c)) for c in np.arange(i, i + 8)]
                ofp.write('\n')
                i += 8

        ofp.close()


