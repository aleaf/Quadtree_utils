'''
example script for writing a usg grid spec file
'''
import os
import sys
sys.path.append('D:ATLData/Documents/GitHub/Quadtree_utils/')

import gridspec as gs

tsf_file = 'D:/ATLData/Wbasin/grid02qtg.tsf'
path = 'D:/ATLData/Wbasin'
topfiles = [os.path.join(path, f) for f in os.listdir(path) if '.top' in f]
botfiles = [os.path.join(path, f) for f in os.listdir(path) if '.bot' in f]
nrow, ncol, nlay = 456, 468, 8
delrc = 5249.3438
xll = -1.12795e+006
yll = 1.71332e+007

# instantiate gridSPC class with parent grid information
spc = gs.gridSPC(nrow, ncol, nlay, xll=xll, yll=yll, dxy=delrc)

# read in top and bottom elevations for each node
spc.read_elevations(topfiles, botfiles)

# read tsf and write grid spec file
spc.write_from_tsf(tsf_file, outfile='D:/ATLData/Wbasin/test.spc')

