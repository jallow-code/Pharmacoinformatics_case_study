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
load results/docking/cross_state_benchmark/WT/best_rmsd_pose.pdb, WT_best
show sticks, WT_best
color orange, WT_best
disable WT_best
load results/docking/cross_state_benchmark/C326Y/best_rmsd_pose.pdb, C326Y_best
show sticks, C326Y_best
color tv_green, C326Y_best
disable C326Y_best
load results/docking/cross_state_benchmark/Y64N/best_rmsd_pose.pdb, Y64N_best
show sticks, Y64N_best
color magenta, Y64N_best
disable Y64N_best
set stick_radius, 0.18, bound_vamifeport
zoom bound_vamifeport, 12
