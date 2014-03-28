import numpy as np
import mfpytools.headfile as hf
import quadtree2structured as qt2s

nodez = np.fromfile('Layer1_nodes.dat', dtype=int, sep=' ')

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
indexfile_basename = 'baseLayer'
hdfile = 'wbasin.out.hds'
headobj = hf.HeadFile(hdfile, gridtype='unstructured')
head = headobj.get_heads_for_all_layers(kstp=1,kper=1)
gridpth = ''
gridname = 'grid02qtg'
nodfile = gridname + '.nod'
refinementlevels = 4

# get number of nodes per layer
nodesperlayfile = gridname + '.nodesperlay.dat'

# Output
indexfile_basename = 'Node_idx' # basename for node index files


parentgrid = qt2s.parentGrid(nrows, ncols, nlay, parent_dxy, units, originX, originY)

index = qt2s.buildIndex(path2indexfiles)

index.QT2base(parentgrid, nodfile)

QTarray = qt2s.QTarray(index)

headobj = hf.HeadFile(hdfile, gridtype='unstructured')
QTheads = headobj.get_heads_for_all_layers(kstp=1, kper=1)

heads = QTarray.mapQT2base(QTheads, parentgrid, base_dxy)

qt2s.Save().array2PDF(heads, 'wbasin_hds.pdf', 'heads', 'ft')