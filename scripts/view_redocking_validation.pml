load data/prepared/proteins/8C03_WT_receptor.pdb, receptor
load data/prepared/ligands/8C03_SZU_bound_reference.pdb, bound_vamifeport
load data/prepared/proteins/8C03_WT_validation.gpf.pdb, docking_box

hide everything
show cartoon, receptor
color slate, receptor
show sticks, bound_vamifeport
color orange, bound_vamifeport
show lines, docking_box
color yellow, docking_box

select pocket_residues, receptor and chain A and resi 43+64+65+68+144+185+314+317+318+320+323+324+325+466+469+470+473+501+504+508
show sticks, pocket_residues
color cyan, pocket_residues

set stick_radius, 0.18, bound_vamifeport
set line_width, 3, docking_box
set cartoon_transparency, 0.35, receptor
zoom bound_vamifeport, 12
