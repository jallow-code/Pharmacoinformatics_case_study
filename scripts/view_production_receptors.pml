load data/prepared/proteins/production/8C03_WT_production.pdb, WT_min
load data/prepared/proteins/production/8C03_C326Y_production.pdb, C326Y_min
load data/prepared/proteins/production/8C03_Y64N_production.pdb, Y64N_min
load data/prepared/proteins/production/8C03_WT_production.gpf.pdb, docking_box
hide everything
show cartoon, WT_min
show cartoon, C326Y_min
show cartoon, Y64N_min
color slate, WT_min
color tv_green, C326Y_min
color magenta, Y64N_min
set cartoon_transparency, 0.35
select pocket, (WT_min or C326Y_min or Y64N_min) and chain A and resi 43+64+65+68+144+185+314+317+318+320+323+324+325+326+466+469+470+473+501+504+508
show sticks, pocket
color cyan, pocket
select mut64, (WT_min or C326Y_min or Y64N_min) and chain A and resi 64
select mut326, (WT_min or C326Y_min or Y64N_min) and chain A and resi 326
show sticks, mut64 or mut326
color yellow, mut64 or mut326
zoom pocket, 12
show lines, docking_box
color orange, docking_box
set line_width, 3, docking_box
