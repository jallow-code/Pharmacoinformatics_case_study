#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
out_dir="$root_dir/data/raw/structures"

mkdir -p "$out_dir"

for pdb_id in 6WBV 8C03 8DL7; do
  url="https://files.rcsb.org/download/${pdb_id}.pdb"
  out_path="$out_dir/${pdb_id}.pdb"
  echo "Downloading ${pdb_id} -> ${out_path}"
  curl -L "$url" -o "$out_path"
done

echo "Done."
