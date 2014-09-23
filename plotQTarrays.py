import os
import numpy as np
import mfpytools.headfile as hf
import quadtree2structured as qt2s


originX = -343800
originY = 5222200
parent_dxy = 1600
base_dxy = 100 # finest discretization in Quadtree grid
nrows = 456
ncols = 468
nlay = 8
units = 'm'

# usg file to plot
hdfile = 'D:/ATLData/Wbasin/pymakemodel/wbasinmodel/wbasin.out.hds'

# grid information
nodfile = 'D:/ATLData/Wbasin/output_usgdata/grid02qtg.nod'
tsfile = 'D:/ATLData/Wbasin/grid02qtg.tsf'
refinementlevels = 4

# index information
indexfile = 'D:/ATLData/Wbasin/pymakemodel/index_files/QTnode_idx.pklz'

# output
outdir = 'D:/ATLData/Wbasin/pymakemodel/output/'

# this class stores input for the parent grid
parentgrid = qt2s.parentGrid(nrows, ncols, nlay, parent_dxy, units, originX, originY)
'''
# make an index file mapping base grid to qt nodes
QTpush = qt2s.qtpush(tsfile, parentgrid, refinementlevels, outfile=indexfile)
QTpush.tsf2inodesg()
'''
# bring in Quadtree head results
headobj = hf.HeadFile(hdfile, gridtype='unstructured')
QTheads = headobj.get_heads_for_all_layers(kstp=1, kper=1)

# instantiate QTarray class
# map heads to base (fine resolution grid), by reading index files
QTarray = qt2s.QTarray(parentgrid, base_dxy, idxfile=indexfile)
base_heads = QTarray.mapQT2base(QTheads, layers=[1])

# save heads to PDF
#QTarray.toPDF('wbasin_hds.pdf', 'heads', 'ft')
outfiles_basename = os.path.join(outdir, 'wbasin_hds')
QTarray.toESRIgrid(outfiles_basename)