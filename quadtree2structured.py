import os
import numpy as np
import mfpytools.headfile as hf
from mfpytools.array_sample import modflow_global_coords 
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages



class parentGrid:

    def __init__(self, nrows, ncols, nlay, parent_dxy, units, originX, originY):
        self.nrows = nrows
        self.ncols = ncols
        self.nlay = nlay
        self.dxy = parent_dxy
        self.units = units
        self.originX = originX
        self.originY = originY


class buildIndex:
    '''
    This class takes information on the parent grid and from the *.nod file, and constructs a "base" uniform grid
    at the finest discretization given in the node file. The base grid is populated with node numbers that correspond
    to the quadtree nodes listed in the nod file. The node number references allow for direct assignment of head values
    to the cells at the base (finest) resolution.
    '''
    def __init__(self, path2indexfiles):
        '''
        Look for pre-existing index files. If an index file for layer 1 is found, making of index files
        will be skipped (saves a lot of time).
        Otherwise the index files will be generated.
        '''
        self.path2indexfiles = path2indexfiles
        self.indexfile_basename = 'baseLayer' # e.g. <basename>1.dat

        if os.path.isfile(os.path.join(path2indexfiles, "{0}1.dat".format(self.indexfile_basename))):
            self.indexfile = True
            print "Index file {0} found, skipping index file creation...".format("{0}1.dat".format(self.indexfile_basename))

        else:
            self.indexfile = False
            print "Index file {0} not found, creating index files...".format("{0}1.dat".format(self.indexfile_basename))


    def QT2base(self, parentGrid, nodfile):

        if not self.indexfile:
            print "reading {0} into pandas dataframe...".format(nodfile)
            self.nodedata = pd.read_csv(nodfile, sep=' ', names=['node', 'lay', 'x', 'y', 'z', 'dx', 'dy', 'dz'])
            self.nodedata['xmin'] = self.nodedata['x'] - self.nodedata['dx'] / 2.
            self.nodedata['ymin'] = self.nodedata['y'] - self.nodedata['dy'] / 2.
            self.nodedata['xmax'] = self.nodedata['x'] + self.nodedata['dx'] / 2.
            self.nodedata['ymax'] = self.nodedata['y'] + self.nodedata['dy'] / 2.
            #self.nodedata['head'] = head

            self.min_spacing = np.min(self.nodedata['dx'])
            print "building base grid at {0} {1} spacing...".format(self.min_spacing, parentGrid.units)
            self.nbaserows = parentGrid.nrows * parentGrid.dxy / self.min_spacing
            self.nbasecols = parentGrid.ncols * parentGrid.dxy / self.min_spacing
            self.X = (np.arange(self.nbasecols) + 0.5) * self.min_spacing + parentGrid.originX
            self.Y = (np.arange(self.nbaserows) + 0.5)[::-1] * self.min_spacing + parentGrid.originY
            uniformgrid = np.empty((self.nbasecols, self.nbaserows))

            for l in range(parentGrid.nlay):
                print "Layer ", l + 1
                # convert dataframe back to numpy to (hopefully) speed iteration
                layernodes = self.nodedata[self.nodedata['lay'] == l+1].to_records()

                #for index, node in layernodes.iterrows():
                knt = 0
                for node in layernodes:
                    knt += 1
                    print "\rUSG node {0}, {1:d}% complete".format(node['node'], 100 * knt / len(layernodes)),
                    basecols = np.where((self.X > node['xmin']) & (self.X < node['xmax']))
                    baserows = np.where((self.Y > node['ymin']) & (self.Y < node['ymax']))
                    #baseinds = np.where((uniformcoords[0] > node['xmin']) & (uniformcoords[0] < node['xmax']) & (uniformcoords[1] > node['ymin']) & (uniformcoords[1] < node['ymax']))
                    baserows = np.squeeze(baserows)
                    baseinds = np.meshgrid(basecols, baserows)
                    uniformgrid[baseinds] = node['node']

                uniformgrid = np.transpose(uniformgrid)
                print "\n\nsaving Layer {1} node references to {0}{1}.dat".format(self.indexfile_basename, l)
                np.savetxt("{0}{1}.dat".format(self.indexfile_basename, l+1), uniformgrid, fmt='%d', delimiter=' ')


class QTarray:

    def __init__(self, buildIndex):

        self.path2indexfiles = buildIndex.path2indexfiles
        if len(self.path2indexfiles) == 0:
            self.path2indexfiles = os.getcwd()

        self.indexfile_basename = buildIndex.indexfile_basename

        self.indexfiles = [f for f in os.listdir(self.path2indexfiles) if self.indexfile_basename in f]

        if len(self.indexfiles) == 0:
            raise IOError('No index files found! Run build_index method first to create index files '
                          'for mapping USG output to regular grid.')

    def mapQT2base(self, QTarray, parentGrid, base_dxy):
        '''
        maps values in MODFLOW-USG array to a base array, using index files generated above
        '''
        nbaserows = parentGrid.nrows * parentGrid.dxy / base_dxy
        nbasecols = parentGrid.ncols * parentGrid.dxy / base_dxy

        print "Mapping Quadtree heads to base grid..."
        uniformgrid = np.empty((parentGrid.nlay, nbaserows, nbasecols))

        for f in self.indexfiles:

            # determine layer index from file name
            layer = int([c for c in f if c.isdigit()][0]) - 1
            print "Layer",
            print " ", layer +1
            # read in index file
            QT_base_indicies = np.fromfile(f, dtype=np.int, sep=' ')

            # set values for current layer from USGarray, base on indicies from index file
            basearray = QTarray[QT_base_indicies]
            uniformgrid[layer, :, :] = np.reshape(basearray, (nbaserows, nbasecols))
            uniformgrid[uniformgrid == 0] = np.nan #convert noflow cells and empty layers to nan

        return uniformgrid


class Save:

    def __init__(self):
        pass

    def array2PDF(self, array, pdfname, title_description=None, zlabel=None, clim=None):

        pdf=PdfPages(pdfname)

        if len(np.shape(array)) == 2:
            nlayers = 1
        else:
            nlayers = np.shape(array)[0]

        print "Saving heads to {0}".format(pdfname)

        for l in range(nlayers):

            # if layer is empty, skip it
            if np.isnan(np.max(array[l, :, :])):
                continue

            print "Layer",
            print " ", l + 1
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_axis_bgcolor('k')
            im = ax.imshow(array[l, :, :], interpolation='nearest')
            if clim:
                im.set_clim(clim)
            cb = fig.colorbar(im)
            cb.set_label(zlabel)
            ax.set_title("Layer {0} {1}".format(l+1, title_description))
            pdf.savefig(fig)
        pdf.close()
'''
class baseArray:

    #This class takes the node index constructed above and assigns values to the base array
s
    def __init__(self, baseNodeIndex):

        if not baseNodeIndex.indexfile:
            nodes =

        hdfile = 'wbasin.out.hds'

nrows = 456
ncols = 468
nlay = 9



baseheads = head[nodez]


'''
'''
Time results on Mac:
Layer 1 with node df converted back to numpy, and searching entire X and Y vector
CPU times: user 2min 23s, sys: 1.87 s, total: 2min 25s
Wall time: 2min 25s
'''
