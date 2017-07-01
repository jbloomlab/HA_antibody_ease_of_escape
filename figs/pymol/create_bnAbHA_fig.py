# Pymol script to create figure of HA with antibodies

# LIST OF PDB FILES
################
# 1RVX: 1934 H1 HA
# 3ZTN: FI6v3 in complex with human H1
# 4GMS: S139/1 in complex with Vic/75 H3
# 4HLZ: C179 in complex with H2
################

import sys

cmd.set('bg_rgb','[1,1,1]') # white
cmd.set('ray_opaque_background','off')
cmd.set('specular', '0')
cmd.set('antialias','2')

cmd.load('1RVX.pdb')
for c in ['G', 'H', 'I', 'J', 'K', 'L']:
	cmd.remove('chain {0}'.format(c))
cmd.select('1RVX_A', '1RVX and chain A')
cmd.select('1RVX_B', '1RVX and chain B')
cmd.select('1RVX_C', '1RVX and chain C')
cmd.color('gray60', '1RVX')

#FI6v3
cmd.load('3ZTN.pdb')
cmd.select('3ZTN_B', '3ZTN and chain B')
cmd.select('3ZTN_HA', '3ZTN and (chain A | chain B)')
cmd.select('FI6v3', '3ZTN and (chain H | chain L)')
cmd.super('3ZTN_B', '1RVX_B')
cmd.color('gray60', '3ZTN_HA')

#C179
cmd.load('4HLZ.pdb')
for c in ['C', 'D', 'E', 'F', 'I', 'J', 'K', 'L']:
	cmd.remove('4HLZ and chain {0}'.format(c))
cmd.select('4HLZ_B', '4HLZ and chain B')
cmd.select('4HLZ_HA', '4HLZ and (chain A | chain B)')
cmd.select('C179', '4HLZ and (chain G | chain H)')
cmd.super('4HLZ_B', '1RVX_B')
cmd.color('gray60', '4HLZ_HA')

#S139/1
cmd.load('4GMS.pdb')
for c in ['C', 'D', 'E', 'F', 'I', 'J', 'M', 'N']:
	cmd.remove('4GMS and chain {0}'.format(c))
cmd.select('4GMS_A', '4GMS and chain A')
cmd.select('4GMS_HA', '4GMS and (chain A | chain B)')
cmd.select('S139', '4GMS and (chain H | chain L)')
cmd.super('4GMS_A', '1RVX_A')
cmd.color('gray60', '4GMS_HA')

cmd.hide('everything')
cmd.show('cartoon')
cmd.show('surface', '1RVX')

cmd.set_view ('\
    -0.067824729,   -0.046455398,   -0.996611416,\
     0.997405410,   -0.027148632,   -0.066612758,\
    -0.023962826,   -0.998546481,    0.048176020,\
     0.000000000,    0.000000000, -616.391479492,\
    69.336593628,   16.585891724,   12.112098694,\
   489.452270508,  743.330932617,  -20.000000000 ')

