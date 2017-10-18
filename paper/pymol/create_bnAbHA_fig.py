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

cmd.load('1RVX_trimer_sequentialnumbering.pdb')
cmd.select('1RVX', '1RVX_trimer_sequentialnumbering')
cmd.select('1RVX_A', '1RVX and chain A')
cmd.select('1RVX_B', '1RVX and chain B')
cmd.select('1RVX_C', '1RVX and chain C')
cmd.color('gray90', '1RVX')

#FI6v3
cmd.load('3ZTN.pdb')
cmd.select('3ZTN_B', '3ZTN and chain B')
cmd.select('3ZTN_HA', '3ZTN and (chain A | chain B)')
cmd.select('FI6v3', '3ZTN and (chain H | chain L)')
cmd.super('3ZTN_B', '1RVX_A')
cmd.remove('3ZTN_HA')

#C179
cmd.load('4HLZ.pdb')
for c in ['C', 'D', 'E', 'F', 'I', 'J', 'K', 'L']:
	cmd.remove('4HLZ and chain {0}'.format(c))
cmd.select('4HLZ_B', '4HLZ and chain B')
cmd.select('4HLZ_HA', '4HLZ and (chain A | chain B)')
cmd.select('C179', '4HLZ and (chain G | chain H)')
cmd.super('4HLZ_B', '1RVX_A')
cmd.remove('4HLZ_HA')

#S139/1
cmd.load('4GMS.pdb')
for c in ['C', 'D', 'E', 'F', 'I', 'J', 'M', 'N']:
	cmd.remove('4GMS and chain {0}'.format(c))
cmd.select('4GMS_A', '4GMS and chain A')
cmd.select('4GMS_HA', '4GMS and (chain A | chain B)')
cmd.select('S139', '4GMS and (chain H | chain L)')
cmd.super('4GMS_A', '1RVX_A')
cmd.remove('4GMS_HA')

cmd.hide('everything')
cmd.show('cartoon')
cmd.show('surface', '1RVX')

cmd.select('HA_Ab', '(1RVX | FI6v3 | S139)')

cmd.color('gold', 'S139')
cmd.color('gadolinium', 'FI6v3')
cmd.color('curium', 'C179')

# H17-L19 footprint
cmd.color('tv_red', '1RVX_A and resi 153')
cmd.color('tv_red', '1RVX_A and resi 154')
cmd.color('tv_red', '1RVX_A and resi 156')
cmd.color('tv_red', '1RVX_A and resi 157')
cmd.color('tv_red', '1RVX_A and resi 158')
cmd.color('tv_red', '1RVX_A and resi 159')
cmd.color('tv_red', '1RVX_A and resi 148')
cmd.color('tv_red', '1RVX_A and resi 136')

# H17-L10 footprint
cmd.color('chocolate', '1RVX and resi 182')
cmd.color('chocolate', '1RVX and resi 186')
cmd.color('chocolate', '1RVX and resi 220')
cmd.color('chocolate', '1RVX and resi 223')
cmd.color('chocolate', '1RVX and resi 234')
cmd.color('chocolate', '1RVX and resi 252')
cmd.color('chocolate', '1RVX and resi 253')
cmd.color('chocolate', '1RVX and resi 254')

# H17-L7 footprint
cmd.color('deepsalmon', '1RVX_B and resi 89')
cmd.color('deepsalmon', '1RVX_B and resi 90')
cmd.color('deepsalmon', '1RVX_B and resi 91')
cmd.color('deepsalmon', '1RVX_B and resi 92')

cmd.set_view ('\
    -0.056619432,   -0.047967870,   -0.997239292,\
     0.997700095,   -0.039890464,   -0.054726347,\
    -0.037155971,   -0.998047054,    0.050115895,\
     0.000000000,   -0.000000000, -557.369567871,\
    79.479003906,    5.570111275,   -3.318444490,\
   430.430236816,  684.308898926,  -20.000000000 ')
