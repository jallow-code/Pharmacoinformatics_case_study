#!/usr/bin/env python
"""Prepare WT re-docking validation assets from the 8C03 complex."""

from __future__ import annotations

import csv
import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW_COMPLEX = ROOT / "data" / "raw" / "structures" / "8C03.pdb"
WT_RECEPTOR_PDB = ROOT / "data" / "prepared" / "proteins" / "8C03_WT_receptor.pdb"

LIGAND_RESNAME = "SZU"
LIGAND_CHAIN = "A"
LIGAND_RESID = 601
PADDING = 4.0
MIN_BOX_SIDE = 18.0
POCKET_CUTOFF = 5.0


@dataclass
class AtomRecord:
    atom_id: int
    atom_name: str
    resname: str
    chain: str
    resid: int
    x: float
    y: float
    z: float
    line: str


def parse_atoms(path: Path) -> list[AtomRecord]:
    atoms: list[AtomRecord] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atoms.append(
                AtomRecord(
                    atom_id=int(line[6:11]),
                    atom_name=line[12:16].strip(),
                    resname=line[17:20].strip(),
                    chain=line[21].strip(),
                    resid=int(line[22:26]),
                    x=float(line[30:38]),
                    y=float(line[38:46]),
                    z=float(line[46:54]),
                    line=line,
                )
            )
    return atoms


def distance(a: AtomRecord, b: AtomRecord) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def extract_ligand_pdb(raw_complex: Path, output_path: Path) -> list[AtomRecord]:
    atoms = parse_atoms(raw_complex)
    ligand_atoms = [
        atom
        for atom in atoms
        if atom.resname == LIGAND_RESNAME
        and atom.chain == LIGAND_CHAIN
        and atom.resid == LIGAND_RESID
    ]
    ligand_ids = {atom.atom_id for atom in ligand_atoms}

    conect_lines: list[str] = []
    with raw_complex.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("CONECT"):
                continue
            fields = line.split()
            atom_ids = [int(x) for x in fields[1:]]
            if atom_ids and atom_ids[0] in ligand_ids and all(x in ligand_ids for x in atom_ids[1:]):
                conect_lines.append(line)

    with output_path.open("w", encoding="utf-8") as handle:
        for atom in ligand_atoms:
            handle.write(atom.line)
        for line in conect_lines:
            handle.write(line)
        handle.write("END\n")

    return ligand_atoms


def compute_box_from_ligand(ligand_atoms: list[AtomRecord]) -> dict[str, float]:
    xs = [atom.x for atom in ligand_atoms]
    ys = [atom.y for atom in ligand_atoms]
    zs = [atom.z for atom in ligand_atoms]

    extents = {
        "x": max(xs) - min(xs),
        "y": max(ys) - min(ys),
        "z": max(zs) - min(zs),
    }
    center = {
        "x": (max(xs) + min(xs)) / 2.0,
        "y": (max(ys) + min(ys)) / 2.0,
        "z": (max(zs) + min(zs)) / 2.0,
    }
    size = {
        axis: max(MIN_BOX_SIDE, extents[axis] + 2.0 * PADDING)
        for axis in ("x", "y", "z")
    }
    return {
        "center_x": round(center["x"], 3),
        "center_y": round(center["y"], 3),
        "center_z": round(center["z"], 3),
        "size_x": round(size["x"], 3),
        "size_y": round(size["y"], 3),
        "size_z": round(size["z"], 3),
        "padding": PADDING,
        "min_box_side": MIN_BOX_SIDE,
    }


def write_box_json(path: Path, params: dict[str, float]) -> None:
    payload = {
        "reference_complex": str(RAW_COMPLEX.relative_to(ROOT)),
        "reference_ligand": f"{LIGAND_RESNAME}:{LIGAND_CHAIN}:{LIGAND_RESID}",
        "box": params,
    }
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")


def write_vina_config(path: Path, receptor_pdbqt: Path, ligand_pdbqt: Path, params: dict[str, float]) -> None:
    text = f"""receptor = {receptor_pdbqt}
ligand = {ligand_pdbqt}

center_x = {params['center_x']}
center_y = {params['center_y']}
center_z = {params['center_z']}

size_x = {params['size_x']}
size_y = {params['size_y']}
size_z = {params['size_z']}

exhaustiveness = 16
num_modes = 20
energy_range = 4
cpu = 0
"""
    path.write_text(text, encoding="utf-8")


def write_pocket_residue_table(path: Path, raw_complex: Path, ligand_atoms: list[AtomRecord]) -> None:
    atoms = parse_atoms(raw_complex)
    protein_atoms = [atom for atom in atoms if atom.line.startswith("ATOM")]

    residue_hits: dict[tuple[str, int, str], float] = {}
    for atom in protein_atoms:
        min_distance = min(distance(atom, lig_atom) for lig_atom in ligand_atoms)
        if min_distance <= POCKET_CUTOFF:
            key = (atom.chain, atom.resid, atom.resname)
            residue_hits[key] = min(min_distance, residue_hits.get(key, min_distance))

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["chain", "resid", "resname", "min_distance_to_bound_ligand"])
        for (chain, resid, resname), min_distance in sorted(residue_hits.items(), key=lambda x: (x[0][0], x[0][1], x[0][2])):
            writer.writerow([chain, resid, resname, f"{min_distance:.3f}"])


def write_pymol_script(path: Path, box_pdb: Path, ligand_pdb: Path, pocket_csv: Path) -> None:
    residues = []
    with pocket_csv.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            residues.append(f"{row['resid']}")
    residue_selection = "+".join(residues)

    text = f"""load {WT_RECEPTOR_PDB.relative_to(ROOT)}, receptor
load {ligand_pdb.relative_to(ROOT)}, bound_vamifeport
load {box_pdb.relative_to(ROOT)}, docking_box

hide everything
show cartoon, receptor
color slate, receptor
show sticks, bound_vamifeport
color orange, bound_vamifeport
show lines, docking_box
color yellow, docking_box

select pocket_residues, receptor and chain A and resi {residue_selection}
show sticks, pocket_residues
color cyan, pocket_residues

set stick_radius, 0.18, bound_vamifeport
set line_width, 3, docking_box
set cartoon_transparency, 0.35, receptor
zoom bound_vamifeport, 12
"""
    path.write_text(text, encoding="utf-8")


def run_command(args: list[str]) -> None:
    subprocess.run(args, check=True, cwd=ROOT)


def main() -> int:
    ligand_dir = ROOT / "data" / "prepared" / "ligands"
    protein_dir = ROOT / "data" / "prepared" / "proteins"
    results_tables = ROOT / "results" / "tables"
    docking_validation_dir = ROOT / "results" / "docking" / "validation"
    config_dir = ROOT / "config"

    ligand_dir.mkdir(parents=True, exist_ok=True)
    protein_dir.mkdir(parents=True, exist_ok=True)
    results_tables.mkdir(parents=True, exist_ok=True)
    docking_validation_dir.mkdir(parents=True, exist_ok=True)
    config_dir.mkdir(parents=True, exist_ok=True)

    ligand_pdb = ligand_dir / "8C03_SZU_bound_reference.pdb"
    ligand_sdf = ligand_dir / "8C03_SZU_bound_reference.sdf"
    ligand_h_sdf = ligand_dir / "8C03_SZU_bound_reference_H.sdf"
    ligand_pdbqt = ligand_dir / "8C03_SZU_bound_reference.pdbqt"

    receptor_prefix = protein_dir / "8C03_WT_validation"
    receptor_pdbqt = protein_dir / "8C03_WT_validation.pdbqt"
    receptor_box_pdb = protein_dir / "8C03_WT_validation.gpf.pdb"

    box_json = config_dir / "redocking_validation_box.json"
    vina_config = config_dir / "vina_redocking.txt"
    pocket_csv = results_tables / "redocking_reference_pocket_residues.csv"
    pymol_script = ROOT / "scripts" / "view_redocking_validation.pml"
    run_script = ROOT / "scripts" / "run_redocking_validation.sh"

    ligand_atoms = extract_ligand_pdb(RAW_COMPLEX, ligand_pdb)
    params = compute_box_from_ligand(ligand_atoms)
    write_box_json(box_json, params)
    write_pocket_residue_table(pocket_csv, RAW_COMPLEX, ligand_atoms)

    run_command(["obabel", str(ligand_pdb), "-O", str(ligand_sdf)])
    run_command(["obabel", str(ligand_sdf), "-O", str(ligand_h_sdf), "-h"])
    run_command(["mk_prepare_ligand.py", "-i", str(ligand_h_sdf), "-o", str(ligand_pdbqt)])
    run_command([
        "mk_prepare_receptor.py",
        "--pdb",
        str(WT_RECEPTOR_PDB),
        "-o",
        str(receptor_prefix),
        "--box_center",
        str(params["center_x"]),
        str(params["center_y"]),
        str(params["center_z"]),
        "--box_size",
        str(params["size_x"]),
        str(params["size_y"]),
        str(params["size_z"]),
    ])

    write_vina_config(vina_config, receptor_pdbqt.relative_to(ROOT), ligand_pdbqt.relative_to(ROOT), params)
    write_pymol_script(pymol_script, receptor_box_pdb, ligand_pdb, pocket_csv)

    run_script.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "${{BASH_SOURCE[0]}}")/.." && pwd)"
cd "$root_dir"

mkdir -p results/docking/validation

vina --config {vina_config.relative_to(ROOT)} \\
  --out results/docking/validation/vamifeport_redock_out.pdbqt \\
  --verbosity 1 \\
  > results/docking/validation/vamifeport_redock.log
""",
        encoding="utf-8",
    )
    run_script.chmod(0o755)

    print("Prepared WT re-docking validation assets:")
    print(f"  ligand pdbqt: {ligand_pdbqt.relative_to(ROOT)}")
    print(f"  receptor pdbqt: {receptor_pdbqt.relative_to(ROOT)}")
    print(f"  box json: {box_json.relative_to(ROOT)}")
    print(f"  vina config: {vina_config.relative_to(ROOT)}")
    print(f"  PyMOL viewer: {pymol_script.relative_to(ROOT)}")
    print(f"  pocket residues: {pocket_csv.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
