import numpy as np
import mfpytools.headfile as hf
import quadtree2structured as qt2s

hdfile = 'wbasin.out.hds'
originX = -343800
originY = 5222200
parent_dxy = 1600
base_dxy = 100 # finest discretization in Quadtree grid
nrows = 456
ncols = 468
nlay = 9
units = 'm'

path2indexfiles = '' # path to index files linking fine-resolution "base" grid with quadtree nodes
indexfile_basename = 'Layer'
hdfile = 'wbasin.out.hds'
gridpth = ''
gridname = 'grid02qtg'
nodfile = gridname + '.nod'
refinementlevels = 4

# get number of nodes per layer
nodesperlayfile = gridname + '.nodesperlay.dat'

# Output
indexfile_basename = 'Node_idx_layer' # basename for node index files

# this class stores input for the parent grid
parentgrid = qt2s.parentGrid(nrows, ncols, nlay, parent_dxy, units, originX, originY)

# instantiate buildIndex class; check for existing index files
index = qt2s.buildIndex(path2indexfiles, indexfile_basename)

# build index arrays at base (finest) resolution if index files don't exist
index.QT2base(parentgrid, nodfile)

# instantiate QTarray class
QTarray = qt2s.QTarray(index)

# bring in Quadtree head results
headobj = hf.HeadFile(hdfile, gridtype='unstructured')
QTheads = headobj.get_heads_for_all_layers(kstp=1, kper=1)

# map heads to base (fine resolution grid), by reading index files
heads = QTarray.mapQT2base(QTheads, parentgrid, base_dxy)

# save heads to PDF
qt2s.Save().array2PDF(heads, 'wbasin_hds.pdf', 'heads', 'ft')
