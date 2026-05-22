#!/usr/bin/env python
"""Dock the validated protonated vamifeport state into WT and mutant production receptors."""

from __future__ import annotations

import csv
import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdFMCS, rdMolAlign


ROOT = Path(__file__).resolve().parents[1]

REFERENCE_LIGAND_PDB = ROOT / "data" / "prepared" / "ligands" / "8C03_SZU_bound_reference.pdb"
REFERENCE_LIGAND_H_SDF = ROOT / "data" / "prepared" / "ligands" / "8C03_SZU_bound_reference_H.sdf"
LIGAND_PDBQT = ROOT / "data" / "prepared" / "ligands" / "redocking_state_panel" / "current_linker_protonated.pdbqt"
PROTEIN_DIR = ROOT / "data" / "prepared" / "proteins" / "production"
BOX_JSON = ROOT / "config" / "redocking_validation_box_optimized.json"

OUTPUT_DIR = ROOT / "results" / "docking" / "cross_state_benchmark"
TABLE_DIR = ROOT / "results" / "tables"
SUMMARY_CSV = TABLE_DIR / "cross_state_benchmark_summary.csv"
POSE_METRICS_CSV = TABLE_DIR / "cross_state_benchmark_pose_metrics.csv"
PYMOL_TOP = ROOT / "scripts" / "view_cross_state_benchmark_top_poses.pml"
PYMOL_BEST = ROOT / "scripts" / "view_cross_state_benchmark_best_rmsd_poses.pml"

CONTACT_CUTOFF = 4.0
NUM_MODES = 20
ENERGY_RANGE = 4
COLORS = {"WT": "orange", "C326Y": "tv_green", "Y64N": "magenta"}

STATE_TO_RECEPTOR = {
    "WT": PROTEIN_DIR / "8C03_WT_production.pdbqt",
    "C326Y": PROTEIN_DIR / "8C03_C326Y_production.pdbqt",
    "Y64N": PROTEIN_DIR / "8C03_Y64N_production.pdbqt",
}

STATE_TO_PDB = {
    "WT": PROTEIN_DIR / "8C03_WT_production.pdb",
    "C326Y": PROTEIN_DIR / "8C03_C326Y_production.pdb",
    "Y64N": PROTEIN_DIR / "8C03_Y64N_production.pdb",
}


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


def load_box_and_vina() -> tuple[dict[str, float], dict[str, int | float]]:
    payload = json.loads(BOX_JSON.read_text(encoding="utf-8"))
    return payload["box"], payload["vina"]


def load_bound_reference() -> tuple[Chem.Mol, Chem.Mol]:
    mol_h = Chem.SDMolSupplier(str(REFERENCE_LIGAND_H_SDF), removeHs=False)[0]
    if mol_h is None:
        raise ValueError(f"Could not load bound reference ligand from {REFERENCE_LIGAND_H_SDF.relative_to(ROOT)}")
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
    if len(reference_match) != reference.GetNumAtoms() or len(pose_match) != pose.GetNumAtoms():
        raise ValueError("Could not build a full heavy-atom mapping for pose RMSD.")
    atom_map = list(zip(pose_match, reference_match))
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


def dock_state(state_id: str, receptor_pdbqt: Path, receptor_pdb: Path, box: dict[str, float], vina_settings: dict[str, int | float], reference_noh: Chem.Mol, reference_contacts: set[tuple[str, int, str]]) -> tuple[list[dict[str, object]], dict[str, object]]:
    receptor_atoms = [atom for atom in parse_atoms(receptor_pdb) if atom.line.startswith("ATOM")]
    result_dir = OUTPUT_DIR / state_id
    result_dir.mkdir(parents=True, exist_ok=True)

    out_pdbqt = result_dir / f"{state_id}_dock_out.pdbqt"
    out_sdf = result_dir / f"{state_id}_dock_out.sdf"
    log_path = result_dir / f"{state_id}_dock.log"

    run_command(
        [
            "vina",
            "--receptor",
            str(receptor_pdbqt),
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
            str(vina_settings["exhaustiveness"]),
            "--num_modes",
            str(vina_settings["num_modes"]),
            "--energy_range",
            str(vina_settings["energy_range"]),
            "--seed",
            str(vina_settings["seed"]),
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
        raise ValueError(f"Pose/score count mismatch for {state_id}: {len(poses)} poses vs {len(scores)} scores")

    pose_rows: list[dict[str, object]] = []
    best_row: dict[str, object] | None = None
    for rank, (pose, score) in enumerate(zip(poses, scores), start=1):
        rmsd = pose_rmsd_against_reference(reference_noh, pose)
        contacts = residue_contacts_from_mol(receptor_atoms, pose, CONTACT_CUTOFF)
        shared_contacts = contacts & reference_contacts
        novel_contacts = contacts - reference_contacts
        mutation_contacts = {entry for entry in contacts if entry[1] in {64, 326}}
        pose_path = result_dir / f"pose_{rank:02d}.pdb"
        write_pose_pdb(pose_path, pose)

        row: dict[str, object] = {
            "state_id": state_id,
            "mode_rank": rank,
            "score_kcal_mol": score,
            "rmsd_to_wt_bound_reference": rmsd,
            "contact_count": len(contacts),
            "shared_contact_count": len(shared_contacts),
            "reference_contact_recovery": len(shared_contacts) / len(reference_contacts) if reference_contacts else 0.0,
            "contacts": format_contacts(contacts),
            "shared_contacts": format_contacts(shared_contacts),
            "novel_contacts": format_contacts(novel_contacts),
            "mutation_contacts": format_contacts(mutation_contacts),
            "pose_pdb": str(pose_path.relative_to(ROOT)),
        }
        pose_rows.append(row)

        if best_row is None or rmsd < float(best_row["rmsd_to_wt_bound_reference"]):
            best_row = row

    top_row = pose_rows[0]
    top_pose_src = ROOT / str(top_row["pose_pdb"])
    best_pose_src = ROOT / str(best_row["pose_pdb"])
    top_pose_dst = result_dir / "top_pose.pdb"
    best_pose_dst = result_dir / "best_rmsd_pose.pdb"
    top_pose_dst.write_text(top_pose_src.read_text(encoding="utf-8"), encoding="utf-8")
    best_pose_dst.write_text(best_pose_src.read_text(encoding="utf-8"), encoding="utf-8")

    summary = {
        "state_id": state_id,
        "top_mode_rank": int(top_row["mode_rank"]),
        "top_score_kcal_mol": float(top_row["score_kcal_mol"]),
        "top_rmsd_to_wt_bound_reference": float(top_row["rmsd_to_wt_bound_reference"]),
        "top_shared_contact_count": int(top_row["shared_contact_count"]),
        "top_reference_contact_recovery": float(top_row["reference_contact_recovery"]),
        "top_shared_contacts": top_row["shared_contacts"],
        "top_mutation_contacts": top_row["mutation_contacts"],
        "top_pose_pdb": str(top_pose_dst.relative_to(ROOT)),
        "best_rmsd_mode_rank": int(best_row["mode_rank"]),
        "best_rmsd_mode_score_kcal_mol": float(best_row["score_kcal_mol"]),
        "best_rmsd_to_wt_bound_reference": float(best_row["rmsd_to_wt_bound_reference"]),
        "best_rmsd_shared_contact_count": int(best_row["shared_contact_count"]),
        "best_rmsd_reference_contact_recovery": float(best_row["reference_contact_recovery"]),
        "best_rmsd_shared_contacts": best_row["shared_contacts"],
        "best_rmsd_mutation_contacts": best_row["mutation_contacts"],
        "best_rmsd_pose_pdb": str(best_pose_dst.relative_to(ROOT)),
        "dock_log": str(log_path.relative_to(ROOT)),
        "dock_sdf": str(out_sdf.relative_to(ROOT)),
    }
    return pose_rows, summary


def write_table(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_pymol_view(path: Path, summaries: list[dict[str, object]], pose_key: str, suffix: str) -> None:
    lines = [
        f"load {STATE_TO_PDB['WT'].relative_to(ROOT)}, WT_receptor",
        f"load {STATE_TO_PDB['C326Y'].relative_to(ROOT)}, C326Y_receptor",
        f"load {STATE_TO_PDB['Y64N'].relative_to(ROOT)}, Y64N_receptor",
        f"load {REFERENCE_LIGAND_PDB.relative_to(ROOT)}, bound_vamifeport",
        "",
        "hide everything",
        "show cartoon, WT_receptor",
        "show cartoon, C326Y_receptor",
        "show cartoon, Y64N_receptor",
        "color slate, WT_receptor",
        "color tv_green, C326Y_receptor",
        "color magenta, Y64N_receptor",
        "show sticks, bound_vamifeport",
        "color orange, bound_vamifeport",
        "set cartoon_transparency, 0.35",
        "select pocket, (WT_receptor or C326Y_receptor or Y64N_receptor) and chain A and resi 43+64+65+68+144+185+314+317+318+320+323+324+325+326+466+469+470+473+501+504+508",
        "show sticks, pocket",
        "color cyan, pocket",
    ]
    for summary in summaries:
        state_id = str(summary["state_id"])
        pose_path = ROOT / str(summary[pose_key])
        object_name = f"{state_id}_{suffix}"
        lines.append(f"load {pose_path.relative_to(ROOT)}, {object_name}")
        lines.append(f"show sticks, {object_name}")
        lines.append(f"color {COLORS[state_id]}, {object_name}")
        lines.append(f"disable {object_name}")
    lines.extend(
        [
            "set stick_radius, 0.18, bound_vamifeport",
            "zoom bound_vamifeport, 12",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    require_file(REFERENCE_LIGAND_PDB)
    require_file(REFERENCE_LIGAND_H_SDF)
    require_file(LIGAND_PDBQT)
    require_file(BOX_JSON)
    for path in list(STATE_TO_RECEPTOR.values()) + list(STATE_TO_PDB.values()):
        require_file(path)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    box, vina_settings = load_box_and_vina()
    reference_ligand_atoms = parse_atoms(REFERENCE_LIGAND_PDB)
    _, reference_noh = load_bound_reference()
    reference_contacts = residue_contacts_from_pdb_atoms(
        [atom for atom in parse_atoms(STATE_TO_PDB["WT"]) if atom.line.startswith("ATOM")],
        reference_ligand_atoms,
        CONTACT_CUTOFF,
    )

    all_pose_rows: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []
    for state_id in ("WT", "C326Y", "Y64N"):
        pose_rows, summary = dock_state(
            state_id,
            STATE_TO_RECEPTOR[state_id],
            STATE_TO_PDB[state_id],
            box,
            vina_settings,
            reference_noh,
            reference_contacts,
        )
        all_pose_rows.extend(pose_rows)
        summaries.append(summary)

    write_table(
        POSE_METRICS_CSV,
        all_pose_rows,
        [
            "state_id",
            "mode_rank",
            "score_kcal_mol",
            "rmsd_to_wt_bound_reference",
            "contact_count",
            "shared_contact_count",
            "reference_contact_recovery",
            "contacts",
            "shared_contacts",
            "novel_contacts",
            "mutation_contacts",
            "pose_pdb",
        ],
    )
    write_table(
        SUMMARY_CSV,
        summaries,
        [
            "state_id",
            "top_mode_rank",
            "top_score_kcal_mol",
            "top_rmsd_to_wt_bound_reference",
            "top_shared_contact_count",
            "top_reference_contact_recovery",
            "top_shared_contacts",
            "top_mutation_contacts",
            "top_pose_pdb",
            "best_rmsd_mode_rank",
            "best_rmsd_mode_score_kcal_mol",
            "best_rmsd_to_wt_bound_reference",
            "best_rmsd_shared_contact_count",
            "best_rmsd_reference_contact_recovery",
            "best_rmsd_shared_contacts",
            "best_rmsd_mutation_contacts",
            "best_rmsd_pose_pdb",
            "dock_log",
            "dock_sdf",
        ],
    )
    write_pymol_view(PYMOL_TOP, summaries, "top_pose_pdb", "top")
    write_pymol_view(PYMOL_BEST, summaries, "best_rmsd_pose_pdb", "best")

    print("Completed cross-state benchmark docking:")
    print(f"  summary: {SUMMARY_CSV.relative_to(ROOT)}")
    print(f"  pose metrics: {POSE_METRICS_CSV.relative_to(ROOT)}")
    print(f"  PyMOL top poses: {PYMOL_TOP.relative_to(ROOT)}")
    print(f"  PyMOL best-rmsd poses: {PYMOL_BEST.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
