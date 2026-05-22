reinitialize
@scripts/presentation_common_settings.pml

load data/prepared/proteins/production/8C03_WT_production.pdb, WT_receptor
load results/docking/analog_panel/vamifeport/WT/best_recovery_pose.pdb, WT_ligand

load data/prepared/proteins/production/8C03_C326Y_production.pdb, C326Y_receptor
load results/docking/analog_panel/vamifeport/C326Y/best_recovery_pose.pdb, C326Y_ligand

load data/prepared/proteins/production/8C03_Y64N_production.pdb, Y64N_receptor
load results/docking/analog_panel/vamifeport/Y64N/best_recovery_pose.pdb, Y64N_ligand

align C326Y_receptor and name CA, WT_receptor and name CA
align C326Y_ligand, WT_ligand
align Y64N_receptor and name CA, WT_receptor and name CA
align Y64N_ligand, WT_ligand

hide everything
show cartoon, WT_receptor or C326Y_receptor or Y64N_receptor
color gray70, WT_receptor or C326Y_receptor or Y64N_receptor

show sticks, WT_ligand or C326Y_ligand or Y64N_ligand
color orange, WT_ligand or C326Y_ligand or Y64N_ligand

select WT_r64, WT_receptor and chain A and resi 64
select WT_r326, WT_receptor and chain A and resi 326
select C326Y_r64, C326Y_receptor and chain A and resi 64
select C326Y_r326, C326Y_receptor and chain A and resi 326
select Y64N_r64, Y64N_receptor and chain A and resi 64
select Y64N_r326, Y64N_receptor and chain A and resi 326

show sticks, WT_r64 or WT_r326 or C326Y_r64 or C326Y_r326 or Y64N_r64 or Y64N_r326
color forest, WT_r64 or C326Y_r64 or Y64N_r64
color magenta, WT_r326 or C326Y_r326 or Y64N_r326

disable C326Y_receptor
disable C326Y_ligand
disable C326Y_r64
disable C326Y_r326
disable Y64N_receptor
disable Y64N_ligand
disable Y64N_r64
disable Y64N_r326

orient WT_ligand or WT_r64 or WT_r326
zoom WT_ligand or WT_r64 or WT_r326, 6
scene WT, store

disable WT_receptor
disable WT_ligand
disable WT_r64
disable WT_r326
enable C326Y_receptor
enable C326Y_ligand
enable C326Y_r64
enable C326Y_r326
scene C326Y, store

disable C326Y_receptor
disable C326Y_ligand
disable C326Y_r64
disable C326Y_r326
enable Y64N_receptor
enable Y64N_ligand
enable Y64N_r64
enable Y64N_r326
scene Y64N, store

scene WT, recall

# Scene switching:
# scene WT, recall
# scene C326Y, recall
# scene Y64N, recall
#
# After choosing a camera angle for one scene:
# get_view
# paste the output into the other scenes with set_view if needed
#
# Render after adjusting the view:
# ray 3200, 2400
# png report/figures/presentation/wt_key_residue_panel.png, dpi=300
