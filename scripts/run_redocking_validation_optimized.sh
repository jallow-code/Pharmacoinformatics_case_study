#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root_dir"

mkdir -p results/docking/validation/optimized

vina --config config/vina_redocking_current_linker_protonated_optimized.txt \
  --out results/docking/validation/optimized/current_linker_protonated_redock_out.pdbqt \
  --verbosity 1 \
  > results/docking/validation/optimized/current_linker_protonated_redock.log
