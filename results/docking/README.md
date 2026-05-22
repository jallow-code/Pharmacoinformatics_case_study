# Docking Results

This directory was intentionally slimmed before repository export.

The retained files are the ones needed to support the report, figure reproduction, and quick structural inspection.

## What was kept

| Subdirectory | Retained contents |
|---|---|
| `analog_panel/` | `best_recovery_pose.pdb` for each ligand-state pair |
| `cross_state_benchmark/` | `top_pose.pdb`, `best_rmsd_pose.pdb`, and docking logs for `WT`, `C326Y`, and `Y64N` |
| `validation/setup_sweep/box_084_ex16/` | final retained validation configuration with representative poses and log |

## What was removed

To keep the repository lighter, the following bulk intermediates were removed:

- per-mode `pose_*.pdb` files from the analog panel
- full `dock_out.pdbqt` and `dock_out.sdf` files from the analog panel
- alternative validation state-panel outputs
- non-selected validation sweep directories
- early one-off WT redocking intermediates that were superseded by the final retained validation setup

The quantitative record of those stages remains in `results/tables/`.
