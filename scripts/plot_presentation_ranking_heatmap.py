#!/usr/bin/env python
"""Create a presentation-ready heatmap for the top-ranked analog panel."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SUMMARY_CSV = ROOT / "results" / "tables" / "analog_cross_state_summary.csv"
RANKED_CSV = ROOT / "results" / "tables" / "final_ranked_hits.csv"
OUT_PNG = ROOT / "report" / "figures" / "presentation" / "ranking_heatmap.png"

STATES = ["WT", "C326Y", "Y64N"]


def load_ranked_top(top_n: int = 6) -> list[str]:
    with RANKED_CSV.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return [row["ligand_id"] for row in rows[:top_n]]


def load_summary() -> dict[tuple[str, str], dict[str, str]]:
    with SUMMARY_CSV.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return {(row["analog_id"], row["state_id"]): row for row in rows}


def main() -> int:
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    top_ids = load_ranked_top(6)
    summary = load_summary()

    recovery = np.zeros((len(top_ids), len(STATES)))
    annotations: list[list[str]] = []

    for i, ligand_id in enumerate(top_ids):
        row_ann = []
        for j, state in enumerate(STATES):
            row = summary[(ligand_id, state)]
            rec = float(row["best_recovery_reference_contact_recovery"])
            rank = int(row["best_recovery_mode_rank"])
            rmsd = float(row["best_recovery_rmsd_to_wt_bound_reference"])
            recovery[i, j] = rec
            row_ann.append(f"{rec:.3f}\nm{rank}\n{rmsd:.2f} Å")
        annotations.append(row_ann)

    fig, ax = plt.subplots(figsize=(8.8, 6.2), dpi=300)
    im = ax.imshow(recovery, cmap="YlGnBu", vmin=0.5, vmax=1.0, aspect="auto")

    ax.set_xticks(range(len(STATES)))
    ax.set_xticklabels(STATES, fontsize=11, fontweight="bold")
    ax.set_yticks(range(len(top_ids)))
    ax.set_yticklabels(top_ids, fontsize=10)
    ax.set_title("Top Analog Panel Across WT, C326Y, and Y64N", fontsize=14, fontweight="bold", pad=12)

    for i in range(len(top_ids)):
        for j in range(len(STATES)):
            value = recovery[i, j]
            color = "black" if value < 0.82 else "white"
            ax.text(
                j,
                i,
                annotations[i][j],
                ha="center",
                va="center",
                fontsize=8.8,
                color=color,
                linespacing=1.15,
            )

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks(np.arange(-0.5, len(STATES), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(top_ids), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cbar.set_label("Reference contact recovery", fontsize=10)
    cbar.ax.tick_params(labelsize=9)

    fig.text(
        0.5,
        0.02,
        "Cell text shows recovery, best-recovery mode rank, and RMSD to the WT bound reference.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    print(f"Wrote {OUT_PNG.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
