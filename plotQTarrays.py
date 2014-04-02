import os
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
tsfile = gridname + '.tsf'
refinementlevels = 4

# get number of nodes per layer
nodesperlayfile = gridname + '.nodesperlay.dat'

# Output
save_node_index = False
indexfile = 'QTnode_idx.pklz'
indexfile_basename = 'Node_idx_layer' # basename for node index files

# this class stores input for the parent grid
parentgrid = qt2s.parentGrid(nrows, ncols, nlay, parent_dxy, units, originX, originY)

if not os.path.isfile(indexfile):
    QTpush = qt2s.qtpush(tsfile, parentgrid, refinementlevels, save_node_index, indexfile)
    inodesg = QTpush.tsf2inodesg()

# bring in Quadtree head results
headobj = hf.HeadFile(hdfile, gridtype='unstructured')
QTheads = headobj.get_heads_for_all_layers(kstp=1, kper=1)

# instantiate QTarray class
# map heads to base (fine resolution grid), by reading index files
QTarray = qt2s.QTarray(index)
heads = QTarray.mapQT2base(QTheads, parentgrid, base_dxy, layers=[1])

# save heads to PDF
qt2s.Save().array2PDF(heads, 'wbasin_hds.pdf', 'heads', 'ft')
'''