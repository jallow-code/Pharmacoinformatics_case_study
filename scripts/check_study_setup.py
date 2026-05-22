#!/usr/bin/env python
"""Validate the repository skeleton for the SLC40A1 study."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "config" / "study.json"


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def load_csv_rows(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    required_paths = [
        CONFIG_PATH,
        ROOT / "data" / "reference" / "selected_variants.csv",
        ROOT / "data" / "reference" / "structures.csv",
        ROOT / "data" / "reference" / "ligands.csv",
        ROOT / "results" / "tables" / "docking_summary.csv",
        ROOT / "results" / "tables" / "md_summary.csv",
        ROOT / "results" / "tables" / "final_ranked_hits.csv",
    ]

    missing = [path for path in required_paths if not path.exists()]
    if missing:
        print("Missing required study files:", file=sys.stderr)
        for path in missing:
            print(f"  - {path}", file=sys.stderr)
        return 1

    config = load_json(CONFIG_PATH)
    states = config.get("states", [])
    structures = load_csv_rows(ROOT / "data" / "reference" / "structures.csv")
    ligands = load_csv_rows(ROOT / "data" / "reference" / "ligands.csv")

    print("Study setup looks structurally complete.")
    print(f"Project: {config.get('project_title', 'unknown')}")
    print(f"States: {', '.join(state['state_id'] for state in states)}")
    print(f"Reference structures: {', '.join(row['pdb_id'] for row in structures)}")
    print(f"Ligand entries: {len(ligands)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
