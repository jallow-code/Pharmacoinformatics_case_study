# SLC40A1 Ferroportin Inhibitor Project

This repository contains a student CADD workflow for `SLC40A1` gain-of-function hemochromatosis.

Core project question:

**Can vamifeport-like small molecules retain plausible inhibitory binding in hepcidin-resistant ferroportin mutants?**

## Start here

If you are opening the project for the first time, read these files in this order:

1. [report/workflow_results_report.qmd](/home/jallow/pharmacoinformatics/case_study/report/workflow_results_report.qmd)
2. [report/presentation_draft.qmd](/home/jallow/pharmacoinformatics/case_study/report/presentation_draft.qmd)
3. [report/README.md](/home/jallow/pharmacoinformatics/case_study/report/README.md)
4. [results/tables/final_ranked_hits.csv](/home/jallow/pharmacoinformatics/case_study/results/tables/final_ranked_hits.csv)
5. [results/tables/analog_cross_state_summary.csv](/home/jallow/pharmacoinformatics/case_study/results/tables/analog_cross_state_summary.csv)

## Study scope

| Item | Choice |
|---|---|
| primary protein states | `WT`, `C326Y`, `Y64N` |
| primary inhibitor modality | small molecules |
| mechanistic comparators | hepcidin and `PR73` |
| structural anchors | `6WBV`, `8C03`, `8DL7` |

## Repository layout

| Directory | Contents |
|---|---|
| `config/` | study manifest, validation box definitions, and docking configuration files |
| `data/reference/` | variant table, structure table, and analog design table |
| `data/raw/` | downloaded starting structures |
| `data/prepared/` | cleaned receptors and ligands used in docking |
| `results/tables/` | main quantitative outputs used in the report |
| `results/docking/` | pose files and docking logs |
| `scripts/` | reproducible workflow scripts and PyMOL scene files |
| `report/` | final report and slide-ready figures |

## Key outputs

| File | Purpose |
|---|---|
| `report/figures/workflow.png` | workflow overview |
| `report/figures/WT_C326Y_Y64N.png` | receptor-state comparison panel |
| `report/figures/validation_pose_overlay.png` | WT docking validation figure |
| `report/figures/analog_modification_pairs.png` | analog design figure |
| `report/figures/top_ranked_analog_panel.png` | top-ranked structures |
| `report/figures/c326y_top_hits_overlay.png` | top `C326Y` structural overlay |
| `results/tables/final_ranked_hits.csv` | final ligand ranking |

## Background note

The original assignment brief is in:

- [task.md](/home/jallow/pharmacoinformatics/case_study/task.md)
