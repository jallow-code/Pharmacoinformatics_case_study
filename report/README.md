# Report Materials

This directory contains the report and the final figures used in it.

## Main files

| File | Purpose |
|---|---|
| `workflow_results_report.qmd` | main technical report |
| `figures/workflow.png` | workflow overview |
| `figures/WT_C326Y_Y64N.png` | `WT`, `C326Y`, and `Y64N` structural panel |
| `figures/validation_pose_overlay.png` | WT validation figure |
| `figures/analog_modification_pairs.png` | analog design panel |
| `figures/top_ranked_analog_panel.png` | top-ranked analog structures |
| `figures/c326y_top_hits_overlay.png` | structural overlay for the top `C326Y` hits |

## Related result tables

The report depends mainly on these tables in `../results/tables/`:

| File | Purpose |
|---|---|
| `final_ranked_hits.csv` | final cross-state ranking |
| `analog_cross_state_summary.csv` | per-state ligand summary |
| `analog_library_manifest.csv` | analog identities and design logic |
| `cross_state_benchmark_summary.csv` | parent benchmark across states |
| `production_receptor_minimization_summary.csv` | restrained minimization summary |
| `redocking_setup_sweep_summary.csv` | WT validation summary |

## Figure source scenes

The PyMOL scene files used to generate the report figures are in `../scripts/`:

- `presentation_mutation_scenes.pml`
- `presentation_validation_pose_overlay.pml`
- `presentation_c326y_top_hits_overlay.pml`
- `presentation_common_settings.pml`
