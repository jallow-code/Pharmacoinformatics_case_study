#!/usr/bin/env python
"""Sweep a small set of WT re-docking parameters for the protonated vamifeport state."""

from __future__ import annotations

import csv
import json
import math
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdFMCS, rdMolAlign


ROOT = Path(__file__).resolve().parents[1]

REFERENCE_LIGAND_PDB = ROOT / "data" / "prepared" / "ligands" / "8C03_SZU_bound_reference.pdb"
REFERENCE_LIGAND_H_SDF = ROOT / "data" / "prepared" / "ligands" / "8C03_SZU_bound_reference_H.sdf"
RECEPTOR_PDB = ROOT / "data" / "prepared" / "proteins" / "8C03_WT_receptor.pdb"
RECEPTOR_PDBQT = ROOT / "data" / "prepared" / "proteins" / "8C03_WT_validation.pdbqt"
BOX_JSON = ROOT / "config" / "redocking_validation_box.json"
BOX_PDB = ROOT / "data" / "prepared" / "proteins" / "8C03_WT_validation.gpf.pdb"
LIGAND_PDBQT = ROOT / "data" / "prepared" / "ligands" / "redocking_state_panel" / "current_linker_protonated.pdbqt"

SWEEP_DIR = ROOT / "results" / "docking" / "validation" / "setup_sweep"
TABLE_DIR = ROOT / "results" / "tables"
SUMMARY_TABLE = TABLE_DIR / "redocking_setup_sweep_summary.csv"
PYMOL_VIEW = ROOT / "scripts" / "view_redocking_setup_sweep_best_configs.pml"

CONTACT_CUTOFF = 4.0
VINA_SEED = 20260521
BOX_SCALES = (1.00, 0.92, 0.84)
EXHAUSTIVENESS_VALUES = (16, 32, 64)
NUM_MODES = 20
ENERGY_RANGE = 4
TOP_CONFIGS_TO_SHOW = 4
COLORS = ["magenta", "marine", "tv_green", "salmon"]


@dataclass(frozen=True)
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


def require_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required input missing: {path.relative_to(ROOT)}")


def run_command(args: list[str], stdout_path: Path | None = None) -> None:
    if stdout_path is None:
        subprocess.run(args, check=True, cwd=ROOT)
        return
    with stdout_path.open("w", encoding="utf-8") as handle:
        subprocess.run(args, check=True, cwd=ROOT, stdout=handle)


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


def distance_xyz(a_xyz: tuple[float, float, float], b_xyz: tuple[float, float, float]) -> float:
    return math.sqrt((a_xyz[0] - b_xyz[0]) ** 2 + (a_xyz[1] - b_xyz[1]) ** 2 + (a_xyz[2] - b_xyz[2]) ** 2)


def load_box() -> dict[str, float]:
    payload = json.loads(BOX_JSON.read_text(encoding="utf-8"))
    return payload["box"]


def load_bound_reference() -> tuple[Chem.Mol, Chem.Mol]:
    mol_h = Chem.SDMolSupplier(str(REFERENCE_LIGAND_H_SDF), removeHs=False)[0]
    if mol_h is None:
        raise ValueError(f"Could not load bound reference from {REFERENCE_LIGAND_H_SDF.relative_to(ROOT)}")
    return mol_h, Chem.RemoveHs(mol_h)


def residue_contacts_from_pdb_atoms(receptor_atoms: list[AtomRecord], ligand_atoms: list[AtomRecord], cutoff: float) -> set[tuple[str, int, str]]:
    hits: set[tuple[str, int, str]] = set()
    ligand_xyz = [(atom.x, atom.y, atom.z) for atom in ligand_atoms if atom.line[76:78].strip() != "H"]
    for atom in receptor_atoms:
        if atom.line[76:78].strip() == "H":
            continue
        atom_xyz = (atom.x, atom.y, atom.z)
        for lig_xyz in ligand_xyz:
            if distance_xyz(atom_xyz, lig_xyz) <= cutoff:
                hits.add((atom.chain, atom.resid, atom.resname))
                break
    return hits


def residue_contacts_from_mol(receptor_atoms: list[AtomRecord], ligand_mol: Chem.Mol, cutoff: float) -> set[tuple[str, int, str]]:
    hits: set[tuple[str, int, str]] = set()
    ligand_noh = Chem.RemoveHs(Chem.Mol(ligand_mol))
    conformer = ligand_noh.GetConformer()
    ligand_xyz = [
        (
            conformer.GetAtomPosition(atom_idx).x,
            conformer.GetAtomPosition(atom_idx).y,
            conformer.GetAtomPosition(atom_idx).z,
        )
        for atom_idx in range(ligand_noh.GetNumAtoms())
    ]
    for atom in receptor_atoms:
        if atom.line[76:78].strip() == "H":
            continue
        atom_xyz = (atom.x, atom.y, atom.z)
        for lig_xyz in ligand_xyz:
            if distance_xyz(atom_xyz, lig_xyz) <= cutoff:
                hits.add((atom.chain, atom.resid, atom.resname))
                break
    return hits


def format_contacts(contacts: set[tuple[str, int, str]]) -> str:
    return ";".join(f"{chain}:{resid}:{resname}" for chain, resid, resname in sorted(contacts))


def pose_rmsd_against_reference(reference_noh: Chem.Mol, pose_mol: Chem.Mol) -> float:
    reference = Chem.Mol(reference_noh)
    pose = Chem.RemoveHs(Chem.Mol(pose_mol))
    mcs = rdFMCS.FindMCS(
        [reference, pose],
        bondCompare=rdFMCS.BondCompare.CompareAny,
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        timeout=10,
    )
    query = Chem.MolFromSmarts(mcs.smartsString)
    reference_match = reference.GetSubstructMatch(query)
    pose_match = pose.GetSubstructMatch(query)
    atom_map = list(zip(pose_match, reference_match))
    if len(atom_map) != reference.GetNumAtoms():
        raise ValueError("Could not map the full heavy-atom graph for RMSD.")
    return rdMolAlign.AlignMol(pose, reference, atomMap=atom_map)


def parse_vina_scores(out_pdbqt: Path) -> list[float]:
    scores: list[float] = []
    with out_pdbqt.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("REMARK VINA RESULT:"):
                scores.append(float(line.split()[3]))
    return scores


def write_pose_pdb(path: Path, mol: Chem.Mol) -> None:
    Chem.MolToPDBFile(mol, str(path))


def scaled_box(base_box: dict[str, float], scale: float) -> dict[str, float]:
    return {
        "center_x": base_box["center_x"],
        "center_y": base_box["center_y"],
        "center_z": base_box["center_z"],
        "size_x": round(base_box["size_x"] * scale, 3),
        "size_y": round(base_box["size_y"] * scale, 3),
        "size_z": round(base_box["size_z"] * scale, 3),
    }


def evaluate_config(config_id: str, box: dict[str, float], exhaustiveness: int, reference_noh: Chem.Mol, receptor_atoms: list[AtomRecord], reference_contacts: set[tuple[str, int, str]]) -> dict[str, object]:
    config_dir = SWEEP_DIR / config_id
    config_dir.mkdir(parents=True, exist_ok=True)

    out_pdbqt = config_dir / f"{config_id}_dock_out.pdbqt"
    out_sdf = config_dir / f"{config_id}_dock_out.sdf"
    log_path = config_dir / f"{config_id}_dock.log"

    run_command(
        [
            "vina",
            "--receptor",
            str(RECEPTOR_PDBQT),
            "--ligand",
            str(LIGAND_PDBQT),
            "--center_x",
            str(box["center_x"]),
            "--center_y",
            str(box["center_y"]),
            "--center_z",
            str(box["center_z"]),
            "--size_x",
            str(box["size_x"]),
            "--size_y",
            str(box["size_y"]),
            "--size_z",
            str(box["size_z"]),
            "--exhaustiveness",
            str(exhaustiveness),
            "--num_modes",
            str(NUM_MODES),
            "--energy_range",
            str(ENERGY_RANGE),
            "--seed",
            str(VINA_SEED),
            "--out",
            str(out_pdbqt),
            "--verbosity",
            "1",
        ],
        stdout_path=log_path,
    )
    run_command(["obabel", str(out_pdbqt), "-O", str(out_sdf)])

    poses = [mol for mol in Chem.SDMolSupplier(str(out_sdf), removeHs=False) if mol is not None]
    scores = parse_vina_scores(out_pdbqt)
    if len(poses) != len(scores):
        raise ValueError(f"Pose/score mismatch for {config_id}: {len(poses)} poses vs {len(scores)} scores")

    top_pose = poses[0]
    top_score = scores[0]
    top_rmsd = pose_rmsd_against_reference(reference_noh, top_pose)
    top_contacts = residue_contacts_from_mol(receptor_atoms, top_pose, CONTACT_CUTOFF)
    top_shared = top_contacts & reference_contacts

    best_rmsd = None
    best_rank = None
    best_score = None
    best_contacts: set[tuple[str, int, str]] | None = None
    best_pose_path = None
    top_pose_path = config_dir / "top_pose.pdb"
    write_pose_pdb(top_pose_path, top_pose)

    for rank, (pose, score) in enumerate(zip(poses, scores), start=1):
        rmsd = pose_rmsd_against_reference(reference_noh, pose)
        if best_rmsd is None or rmsd < best_rmsd:
            best_rmsd = rmsd
            best_rank = rank
            best_score = score
            best_contacts = residue_contacts_from_mol(receptor_atoms, pose, CONTACT_CUTOFF)
            best_pose_path = config_dir / "best_rmsd_pose.pdb"
            write_pose_pdb(best_pose_path, pose)

    assert best_rmsd is not None and best_rank is not None and best_score is not None and best_contacts is not None and best_pose_path is not None
    best_shared = best_contacts & reference_contacts

    return {
        "config_id": config_id,
        "box_scale": box["size_x"] / load_box()["size_x"],
        "size_x": box["size_x"],
        "size_y": box["size_y"],
        "size_z": box["size_z"],
        "exhaustiveness": exhaustiveness,
        "top_score_kcal_mol": top_score,
        "top_rmsd_to_reference": top_rmsd,
        "top_shared_contact_count": len(top_shared),
        "top_reference_contact_recovery": len(top_shared) / len(reference_contacts) if reference_contacts else 0.0,
        "top_shared_contacts": format_contacts(top_shared),
        "best_rmsd_mode_rank": best_rank,
        "best_rmsd_mode_score_kcal_mol": best_score,
        "best_rmsd_to_reference": best_rmsd,
        "best_rmsd_shared_contact_count": len(best_shared),
        "best_rmsd_reference_contact_recovery": len(best_shared) / len(reference_contacts) if reference_contacts else 0.0,
        "best_rmsd_shared_contacts": format_contacts(best_shared),
        "top_pose_pdb": str(top_pose_path.relative_to(ROOT)),
        "best_rmsd_pose_pdb": str(best_pose_path.relative_to(ROOT)),
        "dock_log": str(log_path.relative_to(ROOT)),
        "dock_sdf": str(out_sdf.relative_to(ROOT)),
    }


def write_summary(rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "config_id",
        "box_scale",
        "size_x",
        "size_y",
        "size_z",
        "exhaustiveness",
        "top_score_kcal_mol",
        "top_rmsd_to_reference",
        "top_shared_contact_count",
        "top_reference_contact_recovery",
        "top_shared_contacts",
        "best_rmsd_mode_rank",
        "best_rmsd_mode_score_kcal_mol",
        "best_rmsd_to_reference",
        "best_rmsd_shared_contact_count",
        "best_rmsd_reference_contact_recovery",
        "best_rmsd_shared_contacts",
        "top_pose_pdb",
        "best_rmsd_pose_pdb",
        "dock_log",
        "dock_sdf",
    ]
    with SUMMARY_TABLE.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_pymol_view(rows: list[dict[str, object]]) -> None:
    residue_selection = "43+64+65+68+144+185+314+317+318+320+323+324+325+466+469+470+473+501+504+508"
    lines = [
        f"load {RECEPTOR_PDB.relative_to(ROOT)}, receptor",
        f"load {REFERENCE_LIGAND_PDB.relative_to(ROOT)}, bound_vamifeport",
        f"load {BOX_PDB.relative_to(ROOT)}, docking_box",
        "",
        "hide everything",
        "show cartoon, receptor",
        "color slate, receptor",
        "show sticks, bound_vamifeport",
        "color orange, bound_vamifeport",
        "show lines, docking_box",
        "color yellow, docking_box",
        "",
        f"select pocket_residues, receptor and chain A and resi {residue_selection}",
        "show sticks, pocket_residues",
        "color cyan, pocket_residues",
        "",
    ]
    for color, row in zip(COLORS, rows[:TOP_CONFIGS_TO_SHOW]):
        name = row["config_id"]
        top_pose = ROOT / str(row["top_pose_pdb"])
        best_pose = ROOT / str(row["best_rmsd_pose_pdb"])
        lines.append(f"load {top_pose.relative_to(ROOT)}, {name}_top")
        lines.append(f"show sticks, {name}_top")
        lines.append(f"color {color}, {name}_top")
        lines.append(f"disable {name}_top")
        lines.append(f"load {best_pose.relative_to(ROOT)}, {name}_best")
        lines.append(f"show sticks, {name}_best")
        lines.append(f"color {color}, {name}_best")
        lines.append(f"disable {name}_best")
    lines.extend(
        [
            "",
            "enable bound_vamifeport",
            "set stick_radius, 0.18, bound_vamifeport",
            "set line_width, 3, docking_box",
            "set cartoon_transparency, 0.35, receptor",
            "zoom bound_vamifeport, 12",
        ]
    )
    PYMOL_VIEW.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    for path in [REFERENCE_LIGAND_PDB, REFERENCE_LIGAND_H_SDF, RECEPTOR_PDB, RECEPTOR_PDBQT, BOX_JSON, LIGAND_PDBQT]:
        require_file(path)

    SWEEP_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    base_box = load_box()
    receptor_atoms = [atom for atom in parse_atoms(RECEPTOR_PDB) if atom.line.startswith("ATOM")]
    reference_ligand_atoms = parse_atoms(REFERENCE_LIGAND_PDB)
    reference_contacts = residue_contacts_from_pdb_atoms(receptor_atoms, reference_ligand_atoms, CONTACT_CUTOFF)
    _, reference_noh = load_bound_reference()

    rows: list[dict[str, object]] = []
    for scale in BOX_SCALES:
        box = scaled_box(base_box, scale)
        for exhaustiveness in EXHAUSTIVENESS_VALUES:
            config_id = f"box_{int(round(scale * 100)):03d}_ex{exhaustiveness}"
            rows.append(evaluate_config(config_id, box, exhaustiveness, reference_noh, receptor_atoms, reference_contacts))

    rows.sort(key=lambda row: (float(row["top_rmsd_to_reference"]), int(row["best_rmsd_mode_rank"]), -float(row["top_reference_contact_recovery"])))
    write_summary(rows)
    write_pymol_view(rows)

    print("Completed WT re-docking setup sweep:")
    print(f"  summary: {SUMMARY_TABLE.relative_to(ROOT)}")
    print(f"  PyMOL view: {PYMOL_VIEW.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
