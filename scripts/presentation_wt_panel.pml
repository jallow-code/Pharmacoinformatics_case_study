reinitialize
@scripts/presentation_common_settings.pml

load data/prepared/proteins/production/8C03_WT_production.pdb, receptor
load results/docking/analog_panel/vamifeport/WT/best_recovery_pose.pdb, ligand

hide everything
show cartoon, receptor
color gray70, receptor

select r64, receptor and chain A and resi 64
select r326, receptor and chain A and resi 326
show sticks, ligand or r64 or r326
color orange, ligand
color forest, r64
color magenta, r326

orient ligand or r64 or r326
zoom ligand or r64 or r326, 6

# Adjust the view manually, then render:
# ray 3200, 2400
# png report/figures/presentation/wt_key_residue_panel.png, dpi=300
