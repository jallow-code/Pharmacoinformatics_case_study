#!/usr/bin/env python
"""Build WT and mutant ferroportin receptor models from 8C03."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from openmm.app import PDBFile
from pdbfixer import PDBFixer


DEFAULT_PARENT = Path("data/raw/structures/8C03.pdb")
DEFAULT_OUTPUT_DIR = Path("data/prepared/proteins")
DEFAULT_CHAIN = "A"


def _parse_residue_atoms_from_pdb(path: Path, chain_id: str, residue_id: int) -> dict[str, tuple[float, float, float]]:
    atoms: dict[str, tuple[float, float, float]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            if line[21] != chain_id:
                continue
            if int(line[22:26]) != residue_id:
                continue
            atom_name = line[12:16].strip()
            atoms[atom_name] = (
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            )
    return atoms


def _unit_vector(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    dx = b[0] - a[0]
    dy = b[1] - a[1]
    dz = b[2] - a[2]
    length = math.sqrt(dx * dx + dy * dy + dz * dz)
    return (dx / length, dy / length, dz / length)


def _rewrite_atom_coordinates(line: str, xyz: tuple[float, float, float]) -> str:
    x, y, z = xyz
    return f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"


def repair_y64n_geometry(parent_path: Path, output_path: Path, chain_id: str) -> None:
    """Repair ASN64 terminal geometry using the original TYR64 ring directions.

    PDBFixer correctly changes residue identities and fills missing atoms, but for
    TYR64->ASN the first placed ND2 coordinates can end up too close to CB. This
    function resets OD1 and ND2 using ideal amide bond lengths projected along the
    original TYR CD1/CD2 directions from the parent 8C03 structure.
    """

    parent_atoms = _parse_residue_atoms_from_pdb(parent_path, chain_id, 64)
    model_atoms = _parse_residue_atoms_from_pdb(output_path, chain_id, 64)

    cg_parent = parent_atoms["CG"]
    od1_dir = _unit_vector(cg_parent, parent_atoms["CD1"])
    nd2_dir = _unit_vector(cg_parent, parent_atoms["CD2"])

    cg_model = model_atoms["CG"]
    od1_xyz = (
        cg_model[0] + 1.23 * od1_dir[0],
        cg_model[1] + 1.23 * od1_dir[1],
        cg_model[2] + 1.23 * od1_dir[2],
    )
    nd2_xyz = (
        cg_model[0] + 1.33 * nd2_dir[0],
        cg_model[1] + 1.33 * nd2_dir[1],
        cg_model[2] + 1.33 * nd2_dir[2],
    )

    rewritten_lines: list[str] = []
    with output_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("ATOM") and line[21] == chain_id and int(line[22:26]) == 64:
                atom_name = line[12:16].strip()
                if atom_name == "OD1":
                    line = _rewrite_atom_coordinates(line, od1_xyz)
                elif atom_name == "ND2":
                    line = _rewrite_atom_coordinates(line, nd2_xyz)
            rewritten_lines.append(line)

    with output_path.open("w", encoding="utf-8") as handle:
        handle.writelines(rewritten_lines)


def build_model(
    parent_path: Path,
    output_path: Path,
    chain_id: str,
    mutations: list[str] | None = None,
) -> None:
    fixer = PDBFixer(filename=str(parent_path))
    fixer.removeHeterogens(keepWater=False)

    if mutations:
        fixer.applyMutations(mutations, chain_id)

    # Preserve the experimental backbone and resolved coordinates from 8C03.
    # We only want to complete missing side-chain atoms needed for mutated
    # standard residues, not rebuild the long unresolved loops in the template.
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    with output_path.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build WT and mutant receptor models from the 8C03 template."
    )
    parser.add_argument(
        "--parent",
        type=Path,
        default=DEFAULT_PARENT,
        help="Path to the parent structure PDB file.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for generated receptor models.",
    )
    parser.add_argument(
        "--chain",
        default=DEFAULT_CHAIN,
        help="Protein chain to mutate.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    models = {
        "8C03_WT_receptor.pdb": [],
        "8C03_C326Y_receptor.pdb": ["CYS-326-TYR"],
        "8C03_Y64N_receptor.pdb": ["TYR-64-ASN"],
    }

    for filename, mutations in models.items():
        out_path = args.output_dir / filename
        build_model(args.parent, out_path, args.chain, mutations)
        if "TYR-64-ASN" in mutations:
            repair_y64n_geometry(args.parent, out_path, args.chain)
            print(f"Applied geometry repair to {out_path} for TYR-64-ASN")
        if mutations:
            print(f"Wrote {out_path} with mutations: {', '.join(mutations)}")
        else:
            print(f"Wrote {out_path} (WT receptor)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
