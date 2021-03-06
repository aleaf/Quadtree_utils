import os
import re
import numpy as np
import mfpytools.headfile as hf
from mfpytools.array_sample import modflow_global_coords 
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages



class gridInfo(object):

    def __init__(self, nrows, ncols, nlay, parent_dxy, base_dxy, units, originX, originY, idxfile=None):
        self.nrows = nrows
        self.ncols = ncols
        self.nlay = nlay
        self.dxy = parent_dxy
        self.base_dxy = base_dxy
        self.units = units
        self.originX = originX
        self.originY = originY

        self.nbaserows = self.nrows * self.dxy / self.base_dxy
        self.nbasecols = self.ncols * self.dxy / self.base_dxy

        if not idxfile:
            raise IOError('No index files found! Run build_index method first to create index files '
                          'for mapping USG output to regular grid.')
        else:
            self.idxfile = idxfile


class qtpush:

    def __init__(self, tsfile, parentGrid, levqtmax, bin=True, outfile=None):
        self.tsfile = tsfile
        self.nlaybg = parentGrid.nlay
        self.nrowbg = parentGrid.nrows
        self.ncolbg = parentGrid.ncols
        self.levqtmax = levqtmax
        self.bin = bin
        self.outfile = outfile
        
    def rangefinder(self, qcode, ibgzero, jbgzero, levqtmax, nrowbg, ncolbg):
        
        """
        written by Chris Langevin, USGS Office of Groundwater 4/2/2014
        
        Return the sample grid starting and ending rows and columns for this
        (k,i,j) base grid cell and this quadrant code.
        qcode is the quadrant code
        ibg, jbg is the base grid cell indices (should be zero based)
        levqtmax is the maximum refinement level of the quadtree grid
        nrowbg is the number of rows in the base grid
        ncolbg is the number of columns in the base grid
        """
        
        #calculate size of sample grid using 
        #number of cells per base grid (nsgperbg)
        nsgperbg = 2 ** levqtmax
        nrowsg = nsgperbg * nrowbg
        ncolsg = nsgperbg * ncolbg
        
        #finding starting and ending ranges for base grid cell i,j
        isgmin = ibgzero * nsgperbg
        isgmax = isgmin + nsgperbg
        jsgmin = jbgzero * nsgperbg
        jsgmax = jsgmin + nsgperbg
    
        #process qcode    
        qcode = qcode.strip()
        if len(qcode)>0:
            for q in qcode:
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
                    msg = 'Uknown quadrant code: ' + q + ' for cell ' + str(i,j)
                    raise Exception(msg)
        
        #return with ranges
        return (isgmin,isgmax,jsgmin,jsgmax)
                
    def tsf2inodesg(self):
        """
        written by Chris Langevin, USGS Office of Groundwater 4/2/2014
        Construct inodesg, which is an integer sample grid that has the quadtree
        node number in it.
        levqtmax is the maximum refinement level of the quadtree grid
        nlaybg is the number of layers in the base grid
        nrowbg is the number of rows in the base grid
        ncolbg is the number of columns in the base grid
        """
        tsinfo = open(self.tsfile)
        nsgperbg = 2 ** self.levqtmax
        nrowsg = nsgperbg * self.nrowbg
        ncolsg = nsgperbg * self.ncolbg
        inodesg = np.zeros( (self.nlaybg, nrowsg, ncolsg), dtype=np.int)
        irowidx = np.zeros( (nrowsg * ncolsg), dtype=np.int)
        icolidx = np.zeros( (nrowsg * ncolsg), dtype=np.int)
        nodes = tsinfo.readline()
        for line in tsinfo:
            lnlst = line.split()
            n = int(lnlst[0].replace(',', ''))
            idx = lnlst[1]
            k, i, j = eval(idx)
            qcode=''
            if len(lnlst) > 2:
                qcode = lnlst[2]
            (isgmin, isgmax, jsgmin, jsgmax) = self.rangefinder(qcode, i-1, j-1, self.levqtmax,
                                                        self.nrowbg, self.ncolbg)
            inodesg[k-1, isgmin:isgmax, jsgmin:jsgmax] = n-1
        tsinfo.close()
        if self.bin:
            inodesg.dump(self.outfile)
        return inodesg

        
        
class buildIndex:
    '''
    This class takes information on the parent grid and from the *.nod file, and constructs a "base" uniform grid
    at the finest discretization given in the node file. The base grid is populated with node numbers that correspond
    to the quadtree nodes listed in the nod file. The node number references allow for direct assignment of head values
    to the cells at the base (finest) resolution.
    '''
    def __init__(self, path2indexfiles, indexfile_basename):
        '''
        Look for pre-existing index files. If an index file for layer 1 is found, making of index files
        will be skipped (saves a lot of time).
        Otherwise the index files will be generated.
        '''
        self.path2indexfiles = path2indexfiles
        self.indexfile_basename = indexfile_basename

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
            self.nbaserows = parentGrid.nrows * parentGrid.dxy / self.min_spacing
            self.nbasecols = parentGrid.ncols * parentGrid.dxy / self.min_spacing
            self.X = (np.arange(self.nbasecols) + 0.5) * self.min_spacing + parentGrid.originX
            self.Y = (np.arange(self.nbaserows) + 0.5)[::-1] * self.min_spacing + parentGrid.originY

            print "referencing base grid at {0} {1} spacing...".format(self.min_spacing, parentGrid.units)
            for l in range(parentGrid.nlay):
                print "\nLayer ", l + 1
                # convert dataframe back to numpy to (hopefully) speed iteration
                layernodes = self.nodedata[self.nodedata['lay'] == l+1].to_records()
                lnodes = len(layernodes)
                
                #for index, node in layernodes.iterrows():
                knt = 0
                uniformgrid = np.empty((self.nbasecols, self.nbaserows))
                for node in layernodes:
                    knt += 1
                    print "\rUSG node {0}, {1:d}% complete".format(node['node'], 100 * knt / lnodes),
                    
                    # set y search limits
                    basecols = np.where((self.X > node['xmin']) & (self.X < node['xmax']))
                    baserows = np.where((self.Y > node['ymin']) & (self.Y < node['ymax']))

                    #baseinds = np.where((uniformcoords[0] > node['xmin']) & (uniformcoords[0] < node['xmax']) & (uniformcoords[1] > node['ymin']) & (uniformcoords[1] < node['ymax']))
                    baserows = np.squeeze(baserows)
                    baseinds = np.meshgrid(basecols, baserows)
                    uniformgrid[baseinds] = node['node'] -1 # zero indexing

                uniformgrid = np.transpose(uniformgrid).astype(int)
                print "\nsaving Layer {1} node references to {0}{1}.dat".format(self.indexfile_basename, l+1)    
                #np.savetxt("{0}{1}.dat".format(self.indexfile_basename, l+1), uniformgrid, fmt='%d', delimiter=' ')
                uniformgrid.dump("{0}{1}.dat".format(self.indexfile_basename, l+1))

class QTarray(gridInfo):

    def mapQT2base(self, QTarray, layers='all'):
        '''
        maps values in MODFLOW-USG array to a base array, using index files generated above
        '''

        print "Mapping Quadtree heads to base grid in layers..."

        self.uniformgrid = np.empty((self.nlay, self.nbaserows, self.nbasecols))

        # load single index file for all layers
        self.idx = np.load(self.idxfile)


        if layers == 'all':
            layers = np.arange(1, self.nlay + 1)

        for l in layers:

            print " ", l,
            QT_base_indicies = self.idx[l-1, :, :].astype(int)

            # set values for current layer from USGarray, based on indicies from index file
            basearray = QTarray[QT_base_indicies]
            self.uniformgrid[l-1, :, :] = np.reshape(basearray, (self.nbaserows, self.nbasecols))
            self.uniformgrid[self.uniformgrid == 0] = np.nan #convert noflow cells and empty layers to nan

        return self.uniformgrid


    def toPDF(self, pdfname, title_description=None, zlabel=None, clim=None):

        pdf = PdfPages(pdfname)

        if len(np.shape(self.uniformgrid)) == 2:
            nlayers = 1
        else:
            nlayers = np.shape(self.uniformgrid)[0]

        print "\nSaving heads in layers to {0}...".format(pdfname)

        for l in range(nlayers):

            # if layer is empty, skip it
            if np.isnan(np.max(self.uniformgrid[l, :, :])):
                continue

            print " ", l + 1,
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_axis_bgcolor('k')
            im = ax.imshow(self.uniformgrid[l, :, :], interpolation='nearest')
            if clim:
                im.set_clim(clim)
            cb = fig.colorbar(im)
            cb.set_label(zlabel)
            ax.set_title("Layer {0} {1}".format(l+1, title_description))
            pdf.savefig(fig)
        pdf.close()

    def toESRIgrid(self, basename):

        self.asc_header = 'ncols         {:.0f}\n' \
                          'nrows         {:.0f}\n' \
                          'xllcorner     {}\n' \
                          'yllcorner     {}\n' \
                          'cellsize      {}\n' \
                          'NODATA_value  -9999\n'.format(self.nbasecols, self.nbaserows, self.originX,
                                                         self.originY, self.base_dxy)

        for l in range(self.nlay):

            outfile = basename + '{}.asc'.format(l+1)

            # if layer is empty, skip it
            if np.isnan(np.max(self.uniformgrid[l, :, :])):
                continue

            print "writing {}".format(outfile)
            ofp = open(outfile, 'w')
            ofp.write(self.asc_header)
            np.savetxt(ofp, self.uniformgrid[l, :, :], delimiter=' ', fmt='%.2f')
            ofp.close()