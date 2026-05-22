load data/prepared/proteins/production/8C03_WT_production.pdb, WT_receptor
load data/prepared/proteins/production/8C03_C326Y_production.pdb, C326Y_receptor
load data/prepared/proteins/production/8C03_Y64N_production.pdb, Y64N_receptor
load data/prepared/ligands/8C03_SZU_bound_reference.pdb, bound_vamifeport

hide everything
show cartoon, WT_receptor
show cartoon, C326Y_receptor
show cartoon, Y64N_receptor
color slate, WT_receptor
color tv_green, C326Y_receptor
color magenta, Y64N_receptor
show sticks, bound_vamifeport
color orange, bound_vamifeport
set cartoon_transparency, 0.35
select pocket, (WT_receptor or C326Y_receptor or Y64N_receptor) and chain A and resi 43+64+65+68+144+185+314+317+318+320+323+324+325+326+466+469+470+473+501+504+508
show sticks, pocket
color cyan, pocket
load results/docking/cross_state_benchmark/WT/top_pose.pdb, WT_top
show sticks, WT_top
color orange, WT_top
disable WT_top
load results/docking/cross_state_benchmark/C326Y/top_pose.pdb, C326Y_top
show sticks, C326Y_top
color tv_green, C326Y_top
disable C326Y_top
load results/docking/cross_state_benchmark/Y64N/top_pose.pdb, Y64N_top
show sticks, Y64N_top
color magenta, Y64N_top
disable Y64N_top
set stick_radius, 0.18, bound_vamifeport
zoom bound_vamifeport, 12
