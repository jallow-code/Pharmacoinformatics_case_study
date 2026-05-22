load data/prepared/proteins/production/8C03_C326Y_production.pdb, C326Y_receptor
load data/prepared/ligands/8C03_SZU_bound_reference.pdb, bound_vamifeport
hide everything
show cartoon, C326Y_receptor
color tv_green, C326Y_receptor
show sticks, bound_vamifeport
color orange, bound_vamifeport
select pocket, C326Y_receptor and chain A and resi 43+64+65+68+144+185+314+317+318+320+323+324+325+326+466+469+470+473+501+504+508
show sticks, pocket
color cyan, pocket
set cartoon_transparency, 0.35
load results/docking/analog_panel/T2_benzothiazole/C326Y/best_recovery_pose.pdb, T2_benzothiazole_C326Y_best
show sticks, T2_benzothiazole_C326Y_best
color violet, T2_benzothiazole_C326Y_best
disable T2_benzothiazole_C326Y_best
load results/docking/analog_panel/H3_regiofluoro/C326Y/best_recovery_pose.pdb, H3_regiofluoro_C326Y_best
show sticks, H3_regiofluoro_C326Y_best
color yellow, H3_regiofluoro_C326Y_best
disable H3_regiofluoro_C326Y_best
load results/docking/analog_panel/T1_benzoxazole/C326Y/best_recovery_pose.pdb, T1_benzoxazole_C326Y_best
show sticks, T1_benzoxazole_C326Y_best
color marine, T1_benzoxazole_C326Y_best
disable T1_benzoxazole_C326Y_best
load results/docking/analog_panel/vamifeport/C326Y/best_recovery_pose.pdb, vamifeport_C326Y_best
show sticks, vamifeport_C326Y_best
color orange, vamifeport_C326Y_best
disable vamifeport_C326Y_best
load results/docking/analog_panel/T4_aza_benzimidazole/C326Y/best_recovery_pose.pdb, T4_aza_benzimidazole_C326Y_best
show sticks, T4_aza_benzimidazole_C326Y_best
color limon, T4_aza_benzimidazole_C326Y_best
disable T4_aza_benzimidazole_C326Y_best
set stick_radius, 0.18, bound_vamifeport
zoom bound_vamifeport, 12
