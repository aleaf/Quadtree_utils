'''
example script for writing a usg grid spec file
'''

import sys
sys.path.append('/Users/aleaf/Documents/GitHub/Quadtree_utils/')

import gridspec

nod_file = '/Users/aleaf/Documents/ATLData/Wbasin/output_usgdata/grid02qtg.nod'

spc = gridspec.gridSPC(nod_file)

spc.write_byline(8)