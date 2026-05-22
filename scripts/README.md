# Workflow Scripts

This directory is for the reproducible command-line workflow.

Planned script order:

1. target preparation
2. ligand preparation
3. WT vamifeport re-docking setup
4. WT ligand-state re-docking comparison
5. production docking across WT and mutants
6. static interaction analysis
7. MD system preparation
8. MD analysis and summary tables

Current scripts:

- [check_study_setup.py](/home/jallow/pharmacoinformatics/case_study/scripts/check_study_setup.py)
  validates the study manifest and expected input files.
- [build_analog_library.py](/home/jallow/pharmacoinformatics/case_study/scripts/build_analog_library.py)
  assembles the first-round vamifeport analog set, keeps the protonated linker state fixed, generates 3D conformers, and prepares ligand `pdbqt` files.
- [prepare_redocking_validation.py](/home/jallow/pharmacoinformatics/case_study/scripts/prepare_redocking_validation.py)
  prepares the WT receptor, box, and PyMOL validation scene from `8C03`.
- [run_redocking_state_panel.py](/home/jallow/pharmacoinformatics/case_study/scripts/run_redocking_state_panel.py)
  compares a small set of vamifeport protonation/tautomer states against the same WT re-docking box.
- [optimize_redocking_setup.py](/home/jallow/pharmacoinformatics/case_study/scripts/optimize_redocking_setup.py)
  sweeps a small set of WT re-docking box sizes and exhaustiveness values for the selected protonated ligand state.
- [run_redocking_validation_optimized.sh](/home/jallow/pharmacoinformatics/case_study/scripts/run_redocking_validation_optimized.sh)
  reruns the locked optimized WT validation setup for the selected protonated ligand state.
- [prepare_production_receptors.py](/home/jallow/pharmacoinformatics/case_study/scripts/prepare_production_receptors.py)
  builds the WT, `C326Y`, and `Y64N` production receptor set with restrained local minimization and optimized-box `pdbqt` preparation.
- [run_cross_state_benchmark.py](/home/jallow/pharmacoinformatics/case_study/scripts/run_cross_state_benchmark.py)
  docks the validated protonated vamifeport state into the WT and mutant production receptors and summarizes WT-like pose/contact retention.
- [run_analog_docking_panel.py](/home/jallow/pharmacoinformatics/case_study/scripts/run_analog_docking_panel.py)
  docks the first-round analog library into WT, `C326Y`, and `Y64N`, then summarizes state support and a C326Y-prioritized ranking.

The initial utility script in this skeleton is
[check_study_setup.py](/home/jallow/pharmacoinformatics/case_study/scripts/check_study_setup.py),
which validates the manifest and expected input files.
