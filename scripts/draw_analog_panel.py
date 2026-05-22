#!/usr/bin/env python
"""Draw the parent scaffold and first-round analog modifications for meeting use."""

from __future__ import annotations

import csv
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdDepictor, rdFMCS
from rdkit.Chem.Draw import MolsToGridImage


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_CSV = ROOT / "results" / "tables" / "analog_library_manifest.csv"
RANKED_CSV = ROOT / "results" / "tables" / "final_ranked_hits.csv"
FIG_DIR = ROOT / "report" / "figures"
PAIR_SVG = FIG_DIR / "analog_modification_pairs.svg"
TOP_SVG = FIG_DIR / "top_ranked_analog_panel.svg"
PAIR_PNG = FIG_DIR / "analog_modification_pairs.png"
TOP_PNG = FIG_DIR / "top_ranked_analog_panel.png"

def load_manifest() -> list[dict[str, str]]:
    with MANIFEST_CSV.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def load_ranked_ids(top_n: int = 4) -> list[str]:
    with RANKED_CSV.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return [row["ligand_id"] for row in rows[:top_n]]


def prepare_mol(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Could not parse SMILES: {smiles}")
    rdDepictor.Compute2DCoords(mol)
    return mol


def highlight_relative_to_parent(parent: Chem.Mol, analog: Chem.Mol) -> tuple[list[int], list[int], list[int], list[int]]:
    mcs = rdFMCS.FindMCS(
        [parent, analog],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        timeout=10,
    )
    query = Chem.MolFromSmarts(mcs.smartsString)
    if query is None:
        raise ValueError("Could not build MCS query.")

    parent_match = parent.GetSubstructMatch(query)
    analog_match = analog.GetSubstructMatch(query)
    parent_keep = set(parent_match)
    analog_keep = set(analog_match)

    parent_hilite_atoms = [atom.GetIdx() for atom in parent.GetAtoms() if atom.GetIdx() not in parent_keep]
    analog_hilite_atoms = [atom.GetIdx() for atom in analog.GetAtoms() if atom.GetIdx() not in analog_keep]

    parent_hilite_bonds = [
        bond.GetIdx()
        for bond in parent.GetBonds()
        if bond.GetBeginAtomIdx() in parent_hilite_atoms or bond.GetEndAtomIdx() in parent_hilite_atoms
    ]
    analog_hilite_bonds = [
        bond.GetIdx()
        for bond in analog.GetBonds()
        if bond.GetBeginAtomIdx() in analog_hilite_atoms or bond.GetEndAtomIdx() in analog_hilite_atoms
    ]
    return parent_hilite_atoms, parent_hilite_bonds, analog_hilite_atoms, analog_hilite_bonds


def draw_pairs(rows: list[dict[str, str]]) -> None:
    parent_row = next(row for row in rows if row["analog_id"] == "vamifeport")
    parent = prepare_mol(parent_row["canonical_smiles"])

    mols = []
    legends = []
    highlight_atoms = []
    highlight_bonds = []

    for row in rows:
        if row["analog_id"] == "vamifeport":
            continue
        analog = prepare_mol(row["canonical_smiles"])
        pa, pb, aa, ab = highlight_relative_to_parent(parent, analog)
        mols.extend([Chem.Mol(parent), analog])
        legends.extend(
            [
                f"Parent\nfor {row['analog_id']}",
                f"{row['analog_id']}\n{row['edited_region']}",
            ]
        )
        highlight_atoms.extend([pa, aa])
        highlight_bonds.extend([pb, ab])

    svg = MolsToGridImage(
        mols,
        legends=legends,
        molsPerRow=2,
        subImgSize=(320, 220),
        useSVG=True,
        highlightAtomLists=highlight_atoms,
        highlightBondLists=highlight_bonds,
    )
    PAIR_SVG.write_text(svg, encoding="utf-8")
    png = MolsToGridImage(
        mols,
        legends=legends,
        molsPerRow=2,
        subImgSize=(320, 220),
        useSVG=False,
        highlightAtomLists=highlight_atoms,
        highlightBondLists=highlight_bonds,
    )
    png.save(str(PAIR_PNG))


def draw_top_ranked(rows: list[dict[str, str]], top_ids: list[str]) -> None:
    row_map = {row["analog_id"]: row for row in rows}
    parent = prepare_mol(row_map["vamifeport"]["canonical_smiles"])

    mols = []
    legends = []
    highlight_atoms = []
    highlight_bonds = []

    for analog_id in top_ids:
        row = row_map[analog_id]
        analog = prepare_mol(row["canonical_smiles"])
        if analog_id == "vamifeport":
            pa, pb, aa, ab = [], [], [], []
        else:
            pa, pb, aa, ab = highlight_relative_to_parent(parent, analog)
        mols.append(analog)
        legends.append(f"{analog_id}\n{row['edited_region']}")
        highlight_atoms.append(aa)
        highlight_bonds.append(ab)

    svg = MolsToGridImage(
        mols,
        legends=legends,
        molsPerRow=2,
        subImgSize=(340, 240),
        useSVG=True,
        highlightAtomLists=highlight_atoms,
        highlightBondLists=highlight_bonds,
    )
    TOP_SVG.write_text(svg, encoding="utf-8")
    png = MolsToGridImage(
        mols,
        legends=legends,
        molsPerRow=2,
        subImgSize=(340, 240),
        useSVG=False,
        highlightAtomLists=highlight_atoms,
        highlightBondLists=highlight_bonds,
    )
    png.save(str(TOP_PNG))


def main() -> int:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    rows = load_manifest()
    top_ids = load_ranked_ids(4)
    draw_pairs(rows)
    draw_top_ranked(rows, top_ids)
    print("Wrote analog-structure figures:")
    print(f"  {PAIR_SVG.relative_to(ROOT)}")
    print(f"  {TOP_SVG.relative_to(ROOT)}")
    print(f"  {PAIR_PNG.relative_to(ROOT)}")
    print(f"  {TOP_PNG.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
