reinitialize
@scripts/presentation_common_settings.pml

load data/prepared/proteins/production/8C03_C326Y_production.pdb, receptor
load results/docking/analog_panel/T2_benzothiazole/C326Y/best_recovery_pose.pdb, T2_benzothiazole
load results/docking/analog_panel/H3_regiofluoro/C326Y/best_recovery_pose.pdb, H3_regiofluoro
load results/docking/analog_panel/T1_benzoxazole/C326Y/best_recovery_pose.pdb, T1_benzoxazole
load results/docking/analog_panel/vamifeport/C326Y/best_recovery_pose.pdb, vamifeport

hide everything
show cartoon, receptor
color gray70, receptor

select r64, receptor and chain A and resi 64
select r326, receptor and chain A and resi 326
show sticks, r64 or r326
color forest, r64
color magenta, r326

show sticks, T2_benzothiazole or H3_regiofluoro or T1_benzoxazole or vamifeport
color violet, T2_benzothiazole
color yellow, H3_regiofluoro
color marine, T1_benzoxazole
color orange, vamifeport

orient T2_benzothiazole or H3_regiofluoro or T1_benzoxazole or vamifeport or r64 or r326
zoom T2_benzothiazole or H3_regiofluoro or T1_benzoxazole or vamifeport or r64 or r326, 7

# Adjust the view manually, then render:
# ray 3400, 2600
# png report/figures/presentation/c326y_top_hits_overlay.png, dpi=300
