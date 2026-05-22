#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root_dir"

mkdir -p results/docking/validation

vina --config config/vina_redocking.txt \
  --out results/docking/validation/vamifeport_redock_out.pdbqt \
  --verbosity 1 \
  > results/docking/validation/vamifeport_redock.log
