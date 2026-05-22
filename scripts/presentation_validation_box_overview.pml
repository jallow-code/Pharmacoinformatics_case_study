reinitialize
@scripts/presentation_common_settings.pml

load data/prepared/proteins/8C03_WT_receptor.pdb, receptor
load data/prepared/ligands/8C03_SZU_bound_reference.pdb, bound_vamifeport
load results/docking/validation/setup_sweep/box_084_ex16/top_pose.pdb, redocked_pose
load data/prepared/proteins/8C03_WT_validation.gpf.pdb, docking_box

hide everything
show cartoon, receptor
color gray70, receptor

show sticks, bound_vamifeport
color orange, bound_vamifeport

show sticks, redocked_pose
color magenta, redocked_pose

show lines, docking_box
color marine, docking_box

zoom receptor, 18

# Adjust the view manually, then render:
# ray 3200, 2400
# png report/figures/presentation/validation_box_overview.png, dpi=300
